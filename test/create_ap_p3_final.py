import cobra

import d3flux as d3f
from d3flux.core.flux_layouts import render_model
from d3flux import flux_map

from jinja2 import Template

custom_css = \
"""
{% for item in items %}
text#{{ item }} {
    font-weight: 900;
}
{% endfor %}
text.cofactor {
    fill: #778899;
}
"""

css = Template(custom_css).render(items=['succ_c', 'ac_c', 'etoh_c', 'for_c', 'co2_c', 'lac__D_c'])

model = cobra.io.load_json_model('_all_final_name_formula.json')

# model.add_reaction(ecoli.reactions.ALDD2x)
# d3f.update_cofactors(model, ['nadh_c'])


#model.reactions.EX_glc_e.lower_bound = -1
#model.reactions.EX_xyl_e.lower_bound = -6

import itertools

for obj in itertools.chain(model.reactions, model.metabolites):
    try:
        del obj.notes['map_info']['flux']
        # Initialize empty map_info field in object notes
        if 'map_info' not in obj.notes:
            obj.notes['map_info'] = {}
    except KeyError:
        pass



# model.reactions.PPCK.notes['map_info']['group'] = 2
# model.reactions.MDH.notes['map_info']['group'] = 2
# model.reactions.FUM.notes['map_info']['group'] = 2

# model.reactions.PFL.notes['map_info']['group'] = 'ko'
# model.reactions.ACKr.notes['map_info']['group'] = 'ko'

# Just a single example, but many metabolite names could be better aligned
# model.metabolites.get_by_id('f6p_c').notes['map_info']['align'] = 'center left'

# test
from cobra.flux_analysis import pfba
pfba(model)

html = flux_map(model, figsize=(1280,1024), inactive_alpha=0.5, flux_dict={rxn.id: None for rxn in model.reactions})
with open('_all_final_name_formula.html', 'w') as f:
    f.write('<!DOCTYPE html> <html> <head> \
        <title>Test d3flux page</title> \
        <script src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.2/require.min.js" integrity="sha256-Vjusm6Kh2U7/tb6jBh+MOfxnaf2TWsTph34bMKhC1Qc=" crossorigin="anonymous"></script> \
        <script src="https://code.jquery.com/jquery-3.1.1.min.js" integrity="sha256-hVVnYaiADRTO2PzUGmuLJr8BLUSjGIZsDYGmIJLv2b8=" crossorigin="anonymous"></script> \
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous"> \
        </head> <body>')
    f.write(html.data)
    f.write('</body> </html>')