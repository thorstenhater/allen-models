from allensdk.api.queries.biophysical_api import BiophysicalApi

ids = [497232312, 491766131, 472451419] 
bp = BiophysicalApi()
bp.cache_stimulus = True
for id in ids:
    bp.cache_data(id, working_directory=f'model-{id}-new')
