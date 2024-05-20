import os
import re
import time
import logging
from typing import Union

from dotenv import load_dotenv
import redis
from elasticsearch import Elasticsearch
from retry import retry


load_dotenv()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=os.environ['LOGGING_LEVEL'])


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
                            "gwas_id" if index != 'ieu-b' else "gwas_id.keyword": {
                                "value": str(id)
                            }
                        }
                    },
                    {
                        "range": {
                            "p": {
                                "lte": 0.05
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
            index=index,
            body={
                "query": self._build_es_body_query(id, index)
            }
        )
        return result['count']

    @retry(tries=10, delay=10)
    def run_for_single_dataset(self, gwas_id: str) -> tuple[bool, int]:
        """
        For a single GWAS dataset, count docs, then transform and add to Redis by batch
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
        while len(batch := self.get_es_docs_by_batch(index, id, search_after)) > 0:
            phewas = self._build_redis_zadd_mapping(batch)
            self.add_to_redis(gwas_id, phewas)
            search_after = batch[-1]['sort']
            i += 1
            logging.info('Current batch seq {}, next search after {}'.format(i, search_after))
        return True, n_docs

    def get_es_docs_by_batch(self, index: str, id: Union[str, int], search_after: dict) -> list:
        """
        Fetch Elasticsearch documents by batch
        :param index: GWAS dataset prefix
        :param id: the short GWAS ID
        :param search_after: the search_after parameter for Elasticsearch, obtained from the last hit
        :return: list of Elasticsearch documents
        """
        result = self.es.search(
            index=index,
            size=int(os.environ['BATCH_SIZE']),
            query=self._build_es_body_query(id, index),
            sort=[{
                "chr": {"order": "asc"},
                "position": {"order": "asc"},
                "effect_allele": {"order": "asc"},
                "other_allele": {"order": "asc"}
            }],
            search_after=search_after
        )
        return result['hits']['hits']

    @staticmethod
    def _build_redis_zadd_mapping(docs: list) -> dict:
        """
        Build mappings used for Redis ZADD
        :param docs: list of Elasticsearch documents
        :return: dict of mappings used for ZADD, e.g. {'1:12345:A:T': ('PEHIBaCldykbaP579_By', '0.0470002'), '11:23456:G:C': ('PaCldykbaP57EHIB9_5y', '0.4270002')}
        """
        phewas = {}
        for doc in docs:
            _source = doc['_source']
            phewas['{}:{}:{}:{}'.format(_source['chr'], _source['position'], _source['effect_allele'], _source['other_allele'])] = (doc['_id'], _source['p'])
        return phewas

    def add_to_redis(self, gwas_id: str, phewas: dict) -> None:
        """
        Add members and scores to Redis sorted set
        :param gwas_id: the full GWAS ID
        :param phewas: the mapping string for ZADD
        :return: None
        """
        self.redis.select(int(os.environ['REDIS_DB_PHEWAS']))
        pipe = self.redis.pipeline()
        for key, docid_pval in phewas.items():
            pipe.zadd(key, {gwas_id + ',' + docid_pval[0]: docid_pval[1]})
        results = pipe.execute()
        logging.info('Added to redis: ' + str(results.count(True)))

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


if __name__ == '__main__':
    pi = PhewasIndexing()
    while True:
        if len(tasks := pi.list_pending_tasks_in_redis()) > 0:
            for gwas_id in tasks:
                successful, n_docs = pi.run_for_single_dataset(gwas_id)
                pi.report_task_status_to_redis(gwas_id, successful, n_docs)
        time.sleep(60)
