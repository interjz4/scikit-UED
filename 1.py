#!/usr/bin/env python3
#coding:utf-8

import requests

r = requests.get('https://www.douban.com/') # 豆瓣首页

print(r.status_code)
#print(r.text)