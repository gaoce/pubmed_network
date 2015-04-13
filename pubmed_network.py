#!/usr/bin/env python3
'''
Created on Jan 7, 2014

@author: g
'''


import re
import time
import networkx as nx
import requests as rq
import xml.etree.ElementTree as et


class EutilsBase():
    def __init__(self):
        self._url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        self._para = [('email', 'gaoce@coe.neu.edu'), ('tool', 'biopython')]
        self._delay = 0.333


class ESearch(EutilsBase):
    ''' esearch
    '''
    def __init__(self, kwrds):
        super().__init__()
        self._url_search = self._url + 'esearch.fcgi?'
        self._url_post = self._url + 'epost.fcgi?'
        self._url_summary = self._url + 'esummary.fcgi?'
        self._para = self._para + [('db', 'pubmed'), ('retmax', 10000)]
        self._kwrds = kwrds
        self._articles = []

    def _search_articles(self):
        '''search article on pubmed using kwrds, the results are stored in
        _pmid_list'''
        reg = re.compile(r'(\d\d\d\d)')

        # search in pubmed, get PMID for each article
        para_search = self._para.copy()
        para_search += [('term', self._kwrds)]
        r_search = rq.get(self._url_search, params=para_search)

        root = et.fromstring(r_search.content)
        pmid_tags = root.findall('.//IdList/Id')
        pmid_list = [pmid.text for pmid in pmid_tags]

        # post the resulting PMID, parse webenv and query_key
        para_post = self._para.copy()
        para_post.append(('id', ','.join(pmid_list)))
        r_post = rq.post(self._url_post, data=para_post)

        root = et.fromstring(r_post.content)
        query_key = root.find('QueryKey').text
        webenv = root.find('WebEnv').text

        # get esummary of each pmid
        para_summary = self._para.copy()
        para_summary.extend([('webenv', webenv), ('query_key', query_key)])
        r_summary = rq.get(self._url_summary, params=para_summary)

        root = et.fromstring(r_summary.content)
        doc_elements = root.findall('./DocSum')
        for doc in doc_elements:
            pmid = doc.find('./Id').text
            title = doc.find('./Item[@Name="Title"]').text
            journal = doc.find('./Item[@Name="FullJournalName"]').text
            date = reg.findall(doc.find('./Item[@Name="PubDate"]').text)[0]
            self._articles.append((pmid, {'title':   title,
                                          'journal': journal,
                                          'date':    date}))

    def retrieve(self):
        self._search_articles()
        return self._articles


class Elink(EutilsBase):
    ''' elink
    '''

    def __init__(self, pmid_list, linkname='pubmed_pubmed', buffer_size=400):
        super().__init__()
        self._pmid_list = pmid_list
        self._linkname = linkname
        self._buffer_size = buffer_size
        self._url = self._url + 'elink.fcgi?'

    def _get_links(self):
        para = self._para.copy()
        para.extend([ ('db', 'pubmed'), ('dbfrom', 'pubmed')])
        para.extend([('retmax', 500), ('linkname', self._linkname)])

        edge_list = []
        pmid_list = self._pmid_list

        i = 0
        while i < len(pmid_list):
            pmids = pmid_list[i:i + self._buffer_size]

            # display progress
            print('\r{:.2f}%'.format(i / len(pmid_list) * 100), end='')
            i += self._buffer_size

            paraL = para.copy() # loop parameters
            paraL.extend([('id', pmid) for pmid in pmids])

            # make a post request
            time.sleep(self._delay)
            r = rq.post(self._url, data=paraL)

            # parse the returned xml
            root = et.fromstring(r.content)
            linksets = root.findall('./LinkSet')

            # parse link to get edge
            for element in linksets:
                pmid_from = element.find('./IdList/Id').text
                for link in element.findall('./LinkSetDb/Link/Id'):
                    pmid_to = link.text
                    if pmid_to not in pmid_list or pmid_to == pmid_from:
                        continue
                    edge_list.append((pmid_from, pmid_to))

        print('\r{:.2f}%'.format(100))
        return(edge_list)

    def retrieve(self):
        return self._get_links()


def saveGraph(G, name):
    '''save graph in .gexf and .gpickle format
    '''
    nx.write_gexf(G, name + '.gexf', version='1.2draft')
    nx.write_gpickle(G, name + '.gpickle')


if __name__ == '__main__':
    kwrds = '(("time series") OR "time course") AND "gene expression"'
    es = ESearch(kwrds)
    articles = es.retrieve()
    el = Elink([art[0] for art in articles],'pubmed_pubmed_citedin')
    relations = el.retrieve()

    g = nx.DiGraph()
    g.add_nodes_from(articles)
    g.add_edges_from(relations)
    saveGraph(g, 'new')
