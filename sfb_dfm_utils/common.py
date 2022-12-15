import os

cache_dir=os.path.join(os.path.dirname(__file__),'cache')
os.path.exists(cache_dir) or os.makedirs(cache_dir)
