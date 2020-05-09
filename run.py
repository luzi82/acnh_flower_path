import futsu.json
import futsu.storage
import heapq
import math
import os
import path_algo

S6PASS = 0.9999966
S6ERR = 1-S6PASS
FLOAT_CORRECT = 0.00001

cross = path_algo.cross
cross_verify = path_algo.cross_verify
s_to_v_list_dict = path_algo.s_to_v_list_dict
cross_data_list_to_c_to_p_dict = path_algo.cross_data_list_to_c_to_p_dict

def roll(gene0, gene1, g_to_c_dict, old_gene_set=set()):
  c1 = g_to_c_dict[gene1]
  cross_data_list, chance_sum = cross(gene0, gene1)

  bad_cross_data_list = list(filter(lambda i:i['g'] in old_gene_set, cross_data_list))
  if len(bad_cross_data_list) > 0:
    return None

  cross_data_list = list(filter(lambda i:g_to_c_dict[i['g']]==c1,cross_data_list))
  if len(cross_data_list) == 1:
    cross_data = cross_data_list[0]
    if cross_data['g'] == gene1:
      return {
        'g':cross_data['g'],
        'step_depth':0,
        'step_count':0,
        'add_depth':0,
        'method':'self'
      }
    return {
      'g':cross_data['g'],
      'step_depth':chance_sum/cross_data['chance'],
      'step_count':1,
      'add_depth':chance_sum/cross_data['chance'],
      'method':'cross'
    }
  if len(cross_data_list) == 2:
    cross_data_list0 = list(filter(lambda i:i['g']!=gene1,cross_data_list))
    if len(cross_data_list0) != 1: return None
    cross_data = cross_data_list0[0]
    old_gene_set0 = set(old_gene_set)
    old_gene_set0.add(gene1)
    roll0 = roll(gene0, cross_data['g'], g_to_c_dict, old_gene_set0)
    if roll0 == None: return None
    
    chance_list = map(lambda i:i['chance'], cross_data_list)
    chance_sum0 = sum(chance_list)
    step_depth = chance_sum / chance_sum0
    
    step_count = math.ceil(math.log(S6ERR, 1-(cross_data['chance']/chance_sum0)))
    
    return {
      'g':roll0['g'],
      'step_depth':step_depth,
      'step_count':step_count,
      'add_depth':roll0['add_depth']+step_depth*step_count,
      'method':'roll'
    }
  return None

if __name__ == '__main__':

  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("--flower", type=str)
  parser.add_argument("--myth", type=int, default=1)
  args = parser.parse_args()

  flower_data_path = os.path.join('res','{}.json'.format(args.flower))
  args_myth = args.myth != 0

  flower_data = futsu.json.path_to_data(flower_data_path)
  
  gene_data_list = flower_data['gene_data_list']

  input_txt = 'input/{}.txt'.format(args.flower)
  if futsu.storage.is_blob_exist(input_txt):
    own_gene_set = futsu.fs.file_to_string_list(input_txt)
    own_gene_set = set(filter(lambda i:len(i)>0, own_gene_set))
  else:
    own_gene_set = set()
  
  for gene_data in gene_data_list:
    if 's' not in gene_data: gene_data['s'] = 0
    if 'm' not in gene_data: gene_data['m'] = 0
    gene_data['o'] = 1 if gene_data['g'] in own_gene_set else 0
  
  g_to_c_dict = {
    i['g']: i['c']
    for i in gene_data_list
  }

  gene_done_set = set()

  depth_gene_heap = []
  for gene_data in gene_data_list:
    if gene_data['s'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))
    elif args_myth and gene_data['m'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))
    elif gene_data['o'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))

  print(depth_gene_heap)

  depth_gene_list = []
  gene_to_depth_dict = {}
  gene_to_formula_data_list_dict = {}
  tmp_gene_to_depth_dict = {}

  def add_formul_data(formula_data):
    tmp_depth = tmp_gene_to_depth_dict.get(formula_data['product'],float('inf'))
    if formula_data['total_depth'] - tmp_depth > FLOAT_CORRECT: return
    tmp_gene_to_depth_dict[formula_data['product']] = min(tmp_depth,formula_data['total_depth'])
    heapq.heappush(depth_gene_heap,(formula_data['total_depth'],formula_data['product']))
    if formula_data['product'] not in gene_to_formula_data_list_dict:
      gene_to_formula_data_list_dict[formula_data['product']] = []
    gene_to_formula_data_list_dict[formula_data['product']].append(formula_data)

  for gene_data in gene_data_list:
    if gene_data['s'] > 0:
      add_formul_data({
        'product': gene_data['g'],
        'product.color': g_to_c_dict[gene_data['g']],
        'parent_list': None,
        'step_depth': 0,
        'step_count': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'seed',
      })
    elif args_myth and gene_data['m'] > 0:
      add_formul_data({
        'product': gene_data['g'],
        'product.color': g_to_c_dict[gene_data['g']],
        'parent_list': None,
        'step_depth': 0,
        'step_count': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'myth',
      })
    elif gene_data['o'] > 0:
      add_formul_data({
        'product': gene_data['g'],
        'product.color': g_to_c_dict[gene_data['g']],
        'parent_list': None,
        'step_depth': 0,
        'step_count': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'own',
      })

  while depth_gene_heap:
    depth, gene = heapq.heappop(depth_gene_heap)
    if gene in gene_done_set:
      continue

    gene_done_set.add(gene)
    gene_to_depth_dict[gene] = depth
    depth_gene_list.append((depth, gene))

    for gene_done in gene_done_set:
      gene_chance_list, chance_sum = cross(gene, gene_done)

      c_to_gene_set_dict = {}
      for gene_chance in gene_chance_list:
        g = gene_chance['g']
        c = g_to_c_dict[g]
        if c not in c_to_gene_set_dict:
          c_to_gene_set_dict[c] = set()
        c_to_gene_set_dict[c].add(g)

      for gene_chance in gene_chance_list:
        g = gene_chance['g']
        if g in gene_done_set: continue
        c = g_to_c_dict[g]
        gene_set = c_to_gene_set_dict[c]
        if len(gene_set) > 1: continue
        chance = gene_chance['chance']

        add_depth = chance_sum/chance
        total_depth = add_depth + depth

        add_formul_data({
          'product': g,
          'product.color': g_to_c_dict[g],
          'parent_list': (gene_done, gene),
          'step_depth': add_depth,
          'step_count': 1,
          'add_depth': add_depth,
          'total_depth': total_depth,
          'method': 'cross',
        })
        
      roll0 = roll(gene_done, gene, g_to_c_dict)
      if roll0 and roll0['method'] == 'roll':
        g = roll0['g']
        add_depth = roll0['add_depth']
        total_depth = add_depth + depth

        add_formul_data({
          'product': g,
          'product.color': g_to_c_dict[g],
          'parent_list': (gene_done, gene),
          'step_depth': roll0['step_depth'],
          'step_count': roll0['step_count'],
          'add_depth': add_depth,
          'total_depth': total_depth,
          'method': 'roll',
        })

      roll0 = roll(gene, gene_done, g_to_c_dict)
      if roll0 and roll0['method'] == 'roll':
        g = roll0['g']
        add_depth = roll0['add_depth']
        total_depth = add_depth + depth

        add_formul_data({
          'product': g,
          'product.color': g_to_c_dict[g],
          'parent_list': (gene, gene_done),
          'step_depth': roll0['step_depth'],
          'step_count': roll0['step_count'],
          'add_depth': add_depth,
          'total_depth': total_depth,
          'method': 'roll',
        })

    cross_verify_data_list = []
    for gene_done0 in gene_done_set:
      for gene_done1 in gene_done_set:
        if gene_done1 < gene_done0: continue
        pg0, pg1, vg = gene_done0, gene_done1, gene
        cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict)
        if gene_done0 != gene:
          pg0, pg1, vg = gene, gene_done1, gene_done0
          cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict)
        if gene_done1 != gene:
          pg0, pg1, vg = gene_done0, gene, gene_done1
          cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict)
    for cross_verify_data in cross_verify_data_list:
      add_formul_data({
        'product':      cross_verify_data['product'],
        'product.color': g_to_c_dict[cross_verify_data['product']],
        'parent_list':  cross_verify_data['parent_list'],
        'step_depth':   1,
        'step_count':   cross_verify_data['add_depth'],
        'add_depth':    cross_verify_data['add_depth'],
        'total_depth':  depth+cross_verify_data['add_depth'],
        'method':       'cross_verify',
      })

  print(len(depth_gene_list))
  print(depth_gene_list)
  print(list(sorted(map(lambda i:i[1],depth_gene_list))))

  needed_gene_set = set()

  for depth, gene in reversed(depth_gene_list):
    needed = False
    if '1' not in gene:
      needed = True
    if gene in needed_gene_set:
      needed = True
    if not needed: continue
    formula_data_list = gene_to_formula_data_list_dict[gene]
    formula_data_list = filter(lambda i:i['total_depth'] < depth + FLOAT_CORRECT, formula_data_list)
    formula_data_list = list(formula_data_list)
    for formula_data in formula_data_list:
      if formula_data['total_depth'] > (depth + 0.00001): continue
      #formula_data['product.color'] = g_to_c_dict[formula_data['product']]
      print(formula_data)
      if formula_data['parent_list']:
        for parent_gene in formula_data['parent_list']:
          needed_gene_set.add(parent_gene)
