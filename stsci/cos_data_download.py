from astroquery.mast import Observations
from pathlib import Path
import numpy as np

data_dir = Path('./data/')
data_dir.mkdir(exist_ok=True)
query = Observations.query_criteria(obs_id='LDLM40010')
product_list = Observations.get_product_list(query)

# We want the key where it specifies the raw TIME-TAG both RAWTAG_A and RAWTAG_B which happens to be productSubGroupDescription.
print(product_list.keys())
print(product_list['productSubGroupDescription'])

# We are told to filter only the raw TIME-TAG both RAWTAG_A and RAWTAG_B files and since our product_list is a numpy array, it provides an easy filter mechanicism 
rawtag_product_list = np.where((product_list['productSubGroupDescription'] == 'RAWTAG_A') | (product_list['productSubGroupDescription'] == 'RAWTAG_B'))
print(product_list[rawtag_product_list])

downloads = Observations.download_products(product_list[rawtag_product_list], download_dir=str(data_dir) , extension='fits', mrp_only=False, cache=False)
