
# utilities for working with jupyter notebooks

def first_moment(data):
    return sum([(i + 1) * data[i] for i in range(len(data))])/sum(data)

def nb_load_pre_data():
    import json
    return json.loads(open('pre.spats', 'rb').read())
