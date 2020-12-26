import visflux
import os
import json

# Merge JS dependencies into d3flux.js

# Load the enviroment of visflux
template_path = os.path.dirname(visflux.__file__) + os.path.sep + 'templates'

with open(template_path + os.path.sep + 'include.json') as f:
    data = json.load(f)

include_files = data['local']

merged_js = open(template_path + os.path.sep + 'd3flux.js', 'r+')
# Clear exsisting content
merged_js.truncate(0)

merged_js = open(template_path + os.path.sep +'d3flux.js', mode='w', encoding='utf-8')
for include_file, path in include_files.items():
    temp = open(path, mode='r', encoding='utf-8')
    merged_js.write(temp.read())