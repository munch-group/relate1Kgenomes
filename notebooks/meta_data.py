import pandas as pd
from IPython.display import display

pop_names = {
'ASW': 'African Ancestry in SW U',
'ESN': 'Esan in Nigeria',
'GWD': 'Gambian in Western Division – Mandin',
'LWK': 'Luhya in Webuye, Ken',
'MSL': 'Mende in Sierra Leon',
'YRI': 'Yoruba in Ibadan, Nigeri',
'BEB': 'Bengali in Banglades',
'GIH': 'Gujarati Indians in Houston, Texas, USA',
'PJL': 'Punjabi in Lahore, Pakistan',
'CDX': 'Chinese Dai in Xishuangbanna, China',
'CHB': 'Han Chinese in Beijing, Chin',
'CHS': 'Han Chinese Sout',
'JPT': 'Japanese in Tokyo, Japan',
'KHV': 'Kinh in Ho Chi Minh City, Vietna',
'STU': 'Sri Lankan Tamil in the',
'ITU': 'Indian Telugu in the U.K',
'GBR': 'British From England and Scotlan',
'FIN': 'Finnish in Finla',
'IBS': 'Iberian Populations in Spain',
'TSI': 'Toscani in Itali',
'PUR': 'Puerto Rican in Puerto Rico',
'ACB': 'African Caribbean in Barbado',
'MXL': 'Mexican Ancestry in Los Angeles CA U',
'CLM': 'Colombian in Medellín, Colombia',
'PEL': 'Peruvian in Lima Per',
}
african = ['ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI']
east_asian = ['CDX', 'CHB', 'CHS', 'JPT', 'KHV']
south_asian = ['PJL', 'BEB', 'GIH', 'STU', 'ITU']
european = ['GBR', 'FIN', 'IBS', 'TSI']
central_south_american = ['MXL', 'CLM', 'PEL',]
caribiian = ['PUR', 'ACB']

reg_pop_map = {'Africa': african, 'EastAsian': east_asian, 'SouthAsia': south_asian,
               'Europe': european, 'CentralSouthAmerica': central_south_american,
               'Caribia': caribiian}
pop_reg_map = {pop: reg for reg, pops in reg_pop_map.items() for pop in pops}

populations = list(pop_reg_map.keys())
regions = list(reg_pop_map.keys())
print("Loaded: populations, regions, reg_pop_map, reg_pop_map")
print()
df = pd.DataFrame(dict(label=pop_names.keys(), name=pop_names.values(), region=[pop_reg_map[x] for x in pop_names.keys()]))
print(df.to_string(index=False))