import os
import random
import re
import time
import logging
from collections import defaultdict
from typing import Union
from multiprocessing import Process

from dotenv import load_dotenv
import redis
from elasticsearch import Elasticsearch
from retry import retry


load_dotenv()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=os.environ['LOGGING_LEVEL'])
logging.getLogger('elasticsearch').setLevel(logging.WARNING)


class PhewasIndexing:
    def __init__(self):
        _env = os.environ['ENV']
        self.redis = redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=os.environ['REDIS_DB_TASKS'])
        self.es = Elasticsearch([f"{os.environ['ES_HOST_' + _env]}:{os.environ['ES_PORT']}"])

    def list_pending_tasks_in_redis(self) -> list:
        """
        Fetch list of GWAS IDs from Redis
        :return: list of GWAS IDs
        """
        self.redis.select(int(os.environ['REDIS_DB_TASKS']))
        tasks = self.redis.lrange('pending', 0, -1)
        logging.info('Number of pending tasks: ' + str(len(tasks)))
        return [t.decode('ascii') for t in tasks]

    @staticmethod
    def _split_index_and_id(gwas_id: str) -> tuple[str, str]:
        """
        Split the full GWAS ID into prefix and ID
        :param gwas_id: the full GWAS ID
        :returns: GWAS dataset prefix (index), ID
        """
        index, id = re.match(r'^([\w]+-[\w]+)-([\w]+)', gwas_id).groups()
        return index, id

    @staticmethod
    def _build_es_body_query(id: Union[str, int], index: str) -> dict:
        """
        Build Elasticsearch query string
        :param id: the short GWAS ID
        :param index: GWAS dataset prefix
        :return: query string
        """
        return {
            "bool": {
                "must": [
                    {
                        "term": {
                            "gwas_id": {
                                "value": str(id)
                            }
                        }
                    },
                    {
                        "range": {
                            "p": {
                                "lte": 0.01
                            }
                        }
                    }
                ]
            }
        }

    def count_docs(self, index: str, id: Union[str, int]):
        """
        Count the number of Elasticsearch documents
        :param index: GWAS dataset prefix
        :param id: the short GWAS ID
        :return: number of Elasticsearch documents
        """
        result = self.es.count(
            request_timeout=30,
            index=index,
            body={
                "query": self._build_es_body_query(id, index)
            }
        )
        return result['count']

    def get_es_docs_by_chunk(self, index: str, id: Union[str, int], search_after: dict) -> list:
        """
        Fetch Elasticsearch documents by chunk
        :param index: GWAS dataset prefix, e.g. ebi-a
        :param id: the short GWAS ID
        :param search_after: the search_after parameter for Elasticsearch, obtained from the last hit
        :return: list of Elasticsearch documents
        """
        result = self.es.search(
            request_timeout=90,
            index=index,
            size=int(os.environ['CHUNK_SIZE']),
            query=self._build_es_body_query(id, index),
            sort=[{
                "chr" if index != 'ieu-b' else "chr.keyword": {"order": "asc"},
                "position": {"order": "asc"},
                "effect_allele" if index != 'ieu-b' else "effect_allele.keyword": {"order": "asc"},
                "other_allele" if index != 'ieu-b' else "other_allele.keyword": {"order": "asc"}
            }],
            search_after=search_after
        )
        return result['hits']['hits']

    @staticmethod
    def _build_redis_zadd_mappings(docs: list, index: str) -> (dict, dict):
        """
        Build mappings used for Redis ZADD
        :param docs: list of Elasticsearch documents
        :param index: GWAS dataset prefix, e.g. ebi-a
        :return: two dicts of mapping used for ZADD, e.g. {'1': [('12345', 'A', 'T')], '11': [('23456', 'G', 'C')]} and {'1:12345:A:T': ('PEHIBaCldykbaP579_By', '0.0470002'), '11:23456:G:C': ('PaCldykbaP57EHIB9_5y', '0.4270002')}
        """
        cpalleles = defaultdict(list)
        cpalleles_docids = {}
        for doc in docs:
            _source = doc['_source']
            cpalleles[_source['chr']].append((_source['position'], _source['effect_allele'], _source['other_allele']))
            cpalleles_docids['{}:{}:{}:{}'.format(_source['chr'], _source['position'], _source['effect_allele'], _source['other_allele'])] = (index, doc['_id'], _source['p'])
        return cpalleles, cpalleles_docids

    def add_to_redis(self, cpalleles: dict, cpalleles_docids: dict) -> None:
        """
        Add members and scores to Redis sorted set
        :param cpalleles: the cpalleles mapping dict for ZADD
        :param cpalleles_docids: the cpalleles_docids mapping dict for ZADD
        :return: None
        """
        t = time.time()

        self.redis.select(int(os.environ['REDIS_DB_CPALLELES']))
        pipe = self.redis.pipeline()
        for chr, pos_ea_oa_tuples in cpalleles.items():
            pipe.zadd(chr, {'{}:{}:{}'.format(peo[0], peo[1], peo[2]): peo[0] for peo in pos_ea_oa_tuples})
        results_cpalleles = pipe.execute()
        t2 = time.time()

        self.redis.select(int(os.environ['REDIS_DB_DOCIDS']))
        pipe = self.redis.pipeline()
        for chr_pos_ea_oa, index_docid_pval in cpalleles_docids.items():
            pipe.zadd(chr_pos_ea_oa, {index_docid_pval[0] + ',' + index_docid_pval[1]: index_docid_pval[2]})
        results_cpalleles_docids = pipe.execute()

        logging.debug('Added to redis: {} ({} s), {} ({} s)'.format(str(results_cpalleles.count(True)), str(round(t2 - t, 3)), str(results_cpalleles_docids.count(True)), str(round(time.time() - t2, 3))))

    @retry(tries=10, delay=10)
    def run_for_single_dataset(self, gwas_id: str) -> tuple[bool, int]:
        """
        For a single GWAS dataset, count docs, then transform and add to Redis by chunk
        :param gwas_id: the full GWAS ID
        :return: successful or not, number of docs in Elasticsearch
        """
        try:
            index, id = self._split_index_and_id(gwas_id)
        except AttributeError:
            return False, -1
        n_docs = self.count_docs(index, id)
        logging.info('{} has {} docs'.format(gwas_id, n_docs))
        search_after = None
        i = 0
        while len(chunk := self.get_es_docs_by_chunk(index, id, search_after)) > 0:
            cpalleles, cpalleles_docids = self._build_redis_zadd_mappings(chunk, index)
            self.add_to_redis(cpalleles, cpalleles_docids)
            search_after = chunk[-1]['sort']
            i += 1
            logging.debug('Current chunk seq {}, next search after {}'.format(i, search_after))
        return True, n_docs

    def report_task_status_to_redis(self, gwas_id: str, successful: bool, n_docs: int) -> None:
        """
        Remove a task from 'pending' and add it to 'completed' or 'failed'
        :param gwas_id: the full GWAS ID
        :param successful: whether the task is successful or not
        :param n_docs: number of Elasticsearch documents
        :return: None
        """
        self.redis.select(int(os.environ['REDIS_DB_TASKS']))
        self.redis.lrem('pending', 0, gwas_id)
        if successful:
            self.redis.zadd('completed', {gwas_id: n_docs})
            logging.info('Reported {} as completed with {} docs'.format(gwas_id, n_docs))
        else:
            self.redis.zadd('failed', {gwas_id: int(time.time())})
            logging.info('Reported {} as failed'.format(gwas_id))

    def count_members_in_redis(self) -> (int, int):
        """
        DANGER: THE LUA SCRIPT WILL BLOCK ALL OTHER QUERIES TO THE SERVER. THIS MAY TAKE SEVERAL MINUTES. Run a lua script to count the number of members across all keys in a database
        :return: number of chr:pos:ea:oa combinations, number of Elasticsearch document
        """
        t = time.time()
        self.redis.script_flush()
        lua = """
        redis.call('SELECT', ARGV[1])
        local sum = 0
        local keys = redis.call('KEYS', '*')
        for _, key in ipairs(keys) do
            local card = redis.call('ZCARD', key)
            sum = sum + card
        end
        return sum
        """
        script = self.redis.register_script(lua)
        t2 = time.time()
        logging.info('Script registered as {} in {} s'.format(script.sha, str(round(t2 - t, 3))))

        count_cpalleles = script(args=[int(os.environ['REDIS_DB_CPALLELES'])])
        t3 = time.time()
        logging.info('Count {} cpalleles in {} s'.format(str(count_cpalleles), str(round(t3 - t2, 3))))

        count_cpalleles_docids = script(args=[int(os.environ['REDIS_DB_DOCIDS'])])
        t4 = time.time()
        logging.info('Count {} docs in {} s'.format(str(count_cpalleles_docids), str(round(t4 - t3, 3))))

        return count_cpalleles, count_cpalleles_docids


def single_process(proc_id: int, tasks: list) -> None:
    """
    The target function used by multiprocessing.
    :param proc_id: the process ID
    :param tasks: a list of the full GWAS IDs in this process
    :return: None
    """
    logging.info('Process {} started'.format(proc_id))
    tp = time.time()
    pi = PhewasIndexing()
    for gwas_id in tasks:
        tt = time.time()
        successful, n_docs = pi.run_for_single_dataset(gwas_id)
        pi.report_task_status_to_redis(gwas_id, successful, n_docs)
        logging.info('Task {} completed in {} s'.format(gwas_id, str(round(time.time() - tt, 3))))
    logging.info('Process {} ended in {} s'.format(proc_id, str(round(time.time() - tp, 3))))


if __name__ == '__main__':
    pi = PhewasIndexing()
    n_proc = int(os.environ['N_PROC'])
    while True:
        while len(tasks := pi.list_pending_tasks_in_redis()) > 0:
            random.shuffle(tasks)
            tasks_allocated_to_proc = [tasks[i::n_proc] for i in range(n_proc)]  # if n_proc = 3, tasks of [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] will be divided into [[1, 4, 7, 10], [2, 5, 8], [3, 6, 9]]
            processes = []
            for proc_id in range(n_proc):
                if len(tasks_allocated_to_proc[proc_id]) > 0:
                    proc = Process(target=single_process, args=(proc_id, tasks_allocated_to_proc[proc_id]))
                    proc.start()
                    processes.append(proc)
            for proc in processes:
                proc.join()
        time.sleep(60)
