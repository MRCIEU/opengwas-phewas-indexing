import json
import os
import random
import time

from dotenv import load_dotenv
from elasticsearch import Elasticsearch
import redis
import requests
from retry import retry

load_dotenv()


class Benchmark:
    def __init__(self):
        _env = os.environ['ENV']
        self.redis = redis.Redis(host=os.environ['REDIS_HOST_' + _env], port=os.environ['REDIS_PORT_' + _env], password=os.environ['REDIS_PASS_' + _env], db=os.environ['REDIS_DB_CPALLELES'])
        self.es = Elasticsearch([f"{os.environ['ES_HOST_' + _env]}:{os.environ['ES_PORT']}"])
        self.api_url = os.environ['BENCHMARK_API_URL']
        self.api_key = os.environ['BENCHMARK_API_KEY']
        self.session = requests.Session()

    @retry(tries=10, delay=10)
    def query(self, endpoint: str, variant: list, pval: float):
        t = time.time()
        try:
            r = self.session.post(self.api_url + endpoint, headers={
                'X-TEST-MODE-KEY': self.api_key
            }, data={
                'variant': variant,
                'pval': pval
            })
            assert r.status_code == 200
            return (
                set([(asso['chr'], asso['position'], asso['rsid'], asso['id']) for asso in r.json()]),
                time.time() - t,
                r.headers.get('X-PROCESSING-TIME', '')
            )
        except Exception as e:
            print(e)
            return set(), -1

    @staticmethod
    def _build_es_body_query(chr_pos_combinations: list) -> str:
        queries = []
        body = ''
        for chr_pos in chr_pos_combinations:
            chr_pos = chr_pos.split(':')
            queries.append({})
            queries.append({
                "query": {
                    "bool": {
                        "must": [
                            {
                                "term": {
                                    "CHROM": {
                                        "value": chr_pos[0]
                                    }
                                }
                            },
                            {
                                "term": {
                                    "POS": {
                                        "value": chr_pos[1]
                                    }
                                }
                            },
                            {
                                "term": {
                                    "COMMON": {
                                        "value": "1"
                                    }
                                }
                            }
                        ]
                    }
                }
            })

        for query in queries:
            body += '%s \n' % json.dumps(query)

        return body

    def get_rsid_from_es(self, samples_chrpos):
        response = self.es.msearch(
            request_timeout=180,
            index='snp-base-v0.2',
            body=self._build_es_body_query(samples_chrpos)
        )
        rsids = []
        for r in response['responses']:
            if len(r['hits']['hits']) > 0:
                rsids.append(r['hits']['hits'][0]['_id'])
        return set(rsids)

    def sample_variant(self, n_samples: int, max_radius: int):
        samples_chrpos = set()
        samples_cprange = set()

        for chr in self.redis.scan_iter():
            chr = chr.decode('ascii')
            for pos_alleles in self.redis.zrandmember(chr, n_samples * 5):
                pos_alleles = pos_alleles.decode('ascii')
                pos = int(pos_alleles.split(':', 1)[0])
                samples_chrpos.add('{}:{}'.format(chr, pos))
                samples_cprange.add('{}:{}-{}'.format(chr, pos - random.randint(0, max_radius), pos + random.randint(0, max_radius)))

        samples_rsids = self.get_rsid_from_es(samples_chrpos)

        samples_chrpos = random.sample(list(samples_chrpos), n_samples)
        samples_cprange = random.sample(list(samples_cprange), n_samples)
        samples_rsids = random.sample(list(samples_rsids), n_samples)

        print(len(samples_chrpos), len(samples_cprange), len(samples_rsids))

        result = {
            'rsid': samples_rsids,
            'chrpos': samples_chrpos,
            'cprange': samples_cprange
        }

        with open('samples.json', 'w') as f:
            f.write(json.dumps(result))

        return result

    def run_test(self, pool, n_rounds: int, size_per_variant_type: list[int]):
        with open('results.csv', 'w', buffering=1) as f:
            for i in range(n_rounds):
                for type in pool.keys():
                    for size in size_per_variant_type:
                        print(i, type, str(size))
                        variant = random.sample(pool[type], size)
                        r_phewas = self.query('/phewas', variant, 0.01)
                        r_phewas_fast = self.query('/phewas/fast', variant, 0.01)
                        if r_phewas[1] >= 0 and r_phewas_fast[1] >= 0:
                            if (r_phewas[0] - r_phewas_fast[0]) == set():
                                f.write('{},{},{},{},{}\n'.format(type, str(size), str(round(r_phewas[1], 3)), str(round(r_phewas_fast[1], 3)), r_phewas_fast[2]))
                            else:
                                f.write('ERROR {} {} \n'.format(','.join(variant), str(r_phewas[0] - r_phewas_fast[0])))
                        else:
                            f.write('ERROR {} \n'.format(','.join(variant)))


if __name__ == '__main__':
    b = Benchmark()
    # p = b.query('/phewas', ['rs1280200828'], 0.01)
    # pf = b.query('/phewas/fast', ['rs1280200828'], 0.01)
    # print()
    pool = b.sample_variant(300, 100)
    b.run_test(pool, 10, [1, 2, 4, 8])
