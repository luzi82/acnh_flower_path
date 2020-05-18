import collections
import futsu.json
import futsu.storage
import heapq
import math
import os
import path_algo

from fractions import Fraction

#S6PASS = 0.9999966
#S6ERR = 1-S6PASS
FLOAT_CORRECT = 0.00001

cross = path_algo.cross
cross_verify = path_algo.cross_verify
s_to_v_list_dict = path_algo.s_to_v_list_dict
cross_data_list_to_c_to_p_dict = path_algo.cross_data_list_to_c_to_p_dict
cross_self = path_algo.cross_self
count_one = path_algo.count_one

roll = path_algo.roll
#def roll(*_,**__): return []

#def roll(gene0, gene1, g_to_c_dict, old_gene_set=set()):
#  c1 = g_to_c_dict[gene1]
#  cross_data_list, chance_sum = cross(gene0, gene1)
#
#  bad_cross_data_list = list(filter(lambda i:i['g'] in old_gene_set, cross_data_list))
#  if len(bad_cross_data_list) > 0:
#    return None
#
#  cross_data_list = list(filter(lambda i:g_to_c_dict[i['g']]==c1,cross_data_list))
#  if len(cross_data_list) == 1:
#    cross_data = cross_data_list[0]
#    if cross_data['g'] == gene1:
#      return {
#        'g':cross_data['g'],
#        'step_depth':0,
#        'step_count':0,
#        'add_depth':0,
#        'method':'self'
#      }
#    return {
#      'g':cross_data['g'],
#      'step_depth':chance_sum/cross_data['chance'],
#      'step_count':1,
#      'add_depth':chance_sum/cross_data['chance'],
#      'method':'cross'
#    }
#  if len(cross_data_list) == 2:
#    cross_data_list0 = list(filter(lambda i:i['g']!=gene1,cross_data_list))
#    if len(cross_data_list0) != 1: return None
#    cross_data = cross_data_list0[0]
#    old_gene_set0 = set(old_gene_set)
#    old_gene_set0.add(gene1)
#    roll0 = roll(gene0, cross_data['g'], g_to_c_dict, old_gene_set0)
#    if roll0 == None: return None
#    
#    chance_list = map(lambda i:i['chance'], cross_data_list)
#    chance_sum0 = sum(chance_list)
#    step_depth = chance_sum / chance_sum0
#    
#    step_count = math.ceil(math.log(S6ERR, 1-(cross_data['chance']/chance_sum0)))
#    
#    return {
#      'g':roll0['g'],
#      'step_depth':step_depth,
#      'step_count':step_count,
#      'add_depth':roll0['add_depth']+step_depth*step_count,
#      'method':'roll'
#    }
#  return None

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

  g_to_gene_data_dict = {
    i['g']: i
    for i in gene_data_list
  }
  

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
  path_algo_cache = {}

  gene_done_set = set()

  depth_one_gene_heap = []
  for gene_data in gene_data_list:
    if gene_data['s'] > 0:
      heapq.heappush(depth_one_gene_heap,(0,count_one(gene_data['g']),gene_data['g']))
    elif args_myth and gene_data['m'] > 0:
      heapq.heappush(depth_one_gene_heap,(0,count_one(gene_data['g']),gene_data['g']))
    elif gene_data['o'] > 0:
      heapq.heappush(depth_one_gene_heap,(0,count_one(gene_data['g']),gene_data['g']))

  print(depth_one_gene_heap)

  depth_gene_list = []
  gene_to_depth_dict = {}
  gene_to_formula_data_list_dict = {}
  tmp_gene_to_depth_dict = {}

  def add_formul_data(formula_data):
    tmp_depth = tmp_gene_to_depth_dict.get(formula_data['product'],float('inf'))
    if formula_data['total_depth'] - tmp_depth > FLOAT_CORRECT: return
    tmp_gene_to_depth_dict[formula_data['product']] = min(tmp_depth,formula_data['total_depth'])
    heapq.heappush(depth_one_gene_heap,(formula_data['total_depth'],count_one(formula_data['product']),formula_data['product']))
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
        'begin_depth': 0,
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
        'begin_depth': 0,
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
        'begin_depth': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'own',
      })

  while depth_one_gene_heap:
    depth, _, gene = heapq.heappop(depth_one_gene_heap)
    if gene in gene_done_set:
      continue

    gene_done_set.add(gene)
    gene_to_depth_dict[gene] = depth
    depth_gene_list.append((depth, gene))

    #for gene_done in gene_done_set:
    #  gene_chance_list, chance_sum = cross(gene, gene_done)
    #
    #  c_to_gene_set_dict = {}
    #  for gene_chance in gene_chance_list:
    #    g = gene_chance['g']
    #    c = g_to_c_dict[g]
    #    if c not in c_to_gene_set_dict:
    #      c_to_gene_set_dict[c] = set()
    #    c_to_gene_set_dict[c].add(g)
    #
    #  for gene_chance in gene_chance_list:
    #    g = gene_chance['g']
    #    if g in gene_done_set: continue
    #    c = g_to_c_dict[g]
    #    gene_set = c_to_gene_set_dict[c]
    #    if len(gene_set) > 1: continue
    #    chance = gene_chance['chance']
    #
    #    add_depth = chance_sum/chance
    #    total_depth = add_depth + depth
    #
    #    add_formul_data({
    #      'product': g,
    #      'product.color': g_to_c_dict[g],
    #      'parent_list': (gene_done, gene),
    #      'step_depth': add_depth,
    #      'step_count': 1,
    #      'begin_depth': depth,
    #      'add_depth': add_depth,
    #      'total_depth': total_depth,
    #      'method': 'cross',
    #    })

    for gene_done in gene_done_set:
      cross_self_data_list = cross_self(gene, gene_done, g_to_c_dict, path_algo_cache)
      for cross_self_data in cross_self_data_list:
        add_formul_data({
          'product': cross_self_data['product'],
          'product.color': cross_self_data['product.color'],
          'parent_list': cross_self_data['parent_list'],
          'begin_depth': depth,
          'add_depth': cross_self_data['add_depth'],
          'total_depth': depth+cross_self_data['add_depth'],
          'method': cross_self_data['method'],
        })

    #print('ZECKBYZTBW exit cross_self')

    cross_verify_data_list = []
    for gene_done0 in gene_done_set:
      for gene_done1 in gene_done_set:
        if gene_done1 < gene_done0: continue
        pg0, pg1, vg = gene_done0, gene_done1, gene
        cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
        if gene_done0 != gene:
          pg0, pg1, vg = gene, gene_done1, gene_done0
          cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
        if gene_done1 != gene:
          pg0, pg1, vg = gene_done0, gene, gene_done1
          cross_verify_data_list += cross_verify(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
    for cross_verify_data in cross_verify_data_list:
      #tx = cross_verify_data['parent_list']
      #if tuple(sorted(tx)) == ('001','020','210'):
      #  print('TXWQEJVXCD cross_verify_data={cross_verify_data}'.format(cross_verify_data))
      add_formul_data({
        'product':      cross_verify_data['product'],
        'product.color': g_to_c_dict[cross_verify_data['product']],
        'parent_list':  cross_verify_data['parent_list'],
        'step_depth':   1,
        'step_count':   cross_verify_data['add_depth'],
        'begin_depth':  depth,
        'add_depth':    cross_verify_data['add_depth'],
        'total_depth':  depth+cross_verify_data['add_depth'],
        'method':       'cross_verify',
      })

    roll_data_list = []
    for gene_done0 in gene_done_set:
      for gene_done1 in gene_done_set:
        if gene_done1 < gene_done0: continue
        pg0, pg1, vg = gene_done0, gene_done1, gene
        roll_data_list += roll(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
        if gene_done0 != gene:
          pg0, pg1, vg = gene, gene_done1, gene_done0
          roll_data_list += roll(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
        if gene_done1 != gene:
          pg0, pg1, vg = gene_done0, gene, gene_done1
          roll_data_list += roll(pg0, pg1, vg, g_to_c_dict, depth, tmp_gene_to_depth_dict, path_algo_cache)
    for roll_data in roll_data_list:
      add_formul_data({
        'product':      roll_data['product'],
        'product.color': g_to_c_dict[roll_data['product']],
        'parent_list':  roll_data['parent_list'],
        'step_depth':   1,
        'step_count':   roll_data['add_depth'],
        'begin_depth':  depth,
        'add_depth':    roll_data['add_depth'],
        'total_depth':  depth+roll_data['add_depth'],
        'method':       'roll',
      })

  #print(len(depth_gene_list))
  #print(depth_gene_list)
  #print(list(sorted(map(lambda i:i[1],depth_gene_list))))

  good_gene_list = depth_gene_list
  good_gene_list = map(lambda i:i[1], good_gene_list)
  good_gene_list = filter(lambda i:'1' not in i, good_gene_list)
  good_gene_list = sorted(good_gene_list)
  print('QJBVOIFLVS {}'.format(len(good_gene_list)))
  print('MYYEZVMACM {}'.format(good_gene_list))

  needed_gene_set = set()

  s0_formula_data_list = []

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
      print('XECAGVSANZ formula_data={formula_data}'.format(
        formula_data=formula_data
      ))
      s0_formula_data_list.append(formula_data)
      if formula_data['parent_list']:
        for parent_gene in formula_data['parent_list']:
          needed_gene_set.add(parent_gene)

  pg_to_formula_data_list_dict = {}
  for formula_data in s0_formula_data_list:
    pg_to_formula_data_list_dict[formula_data['product']] = []
  for formula_data in s0_formula_data_list:
    pg_to_formula_data_list_dict[formula_data['product']].append(formula_data)

  must_g_set = set()
  must_g_queue = collections.deque()

  # mark 0/2 gene
  for pg in pg_to_formula_data_list_dict:
    if '1' not in pg:
      must_g_set.add(pg)
      must_g_queue.append(pg)

  while len(must_g_queue)>0:
    g = must_g_queue.popleft()
    formula_data_list = pg_to_formula_data_list_dict[g]
    
    parent_list = formula_data_list[0]['parent_list']
    if not parent_list: continue
    common_parent_set = set(parent_list)
    for formula_data in formula_data_list:
      parent_list = formula_data['parent_list']
      if not parent_list:
        common_parent_set = set()
      else:
        common_parent_set = common_parent_set & set(formula_data['parent_list'])
    
    for g0 in common_parent_set:
      if g0 in must_g_set: continue
      must_g_set.add(g0)
      must_g_queue.append(g0)

  print('BQHWSRIMFK must_g_set={must_g_set}'.format(
    must_g_set=must_g_set
  ))

  pg_to_extra_cost_tuple_dict = {}
  for formula_data in reversed(s0_formula_data_list):
    pg = formula_data['product']

    if formula_data['total_depth'] == 0:
      extra_cost_tuple = (0,)
    else:
      extra_cost_tuple = map(lambda i:math.ceil(pg_to_extra_cost_tuple_dict[i][0]), formula_data['parent_list'])
      extra_cost_tuple = list(extra_cost_tuple)
      extra_cost_tuple.append(math.ceil(max(0,min(extra_cost_tuple))+formula_data['add_depth']))
      extra_cost_tuple = tuple(reversed(sorted(extra_cost_tuple)))
    formula_data['extra_cost_tuple'] = extra_cost_tuple

    if g_to_gene_data_dict[pg]['s']>0:
      extra_cost_tuple = (-3,)
    elif args_myth and g_to_gene_data_dict[pg]['m']>0:
      extra_cost_tuple = (-2,)
    elif g_to_gene_data_dict[pg]['o']>0:
      extra_cost_tuple = (-1,)
    elif pg in must_g_set:
      extra_cost_tuple = (0,)
    old_extra_cost_tuple = pg_to_extra_cost_tuple_dict.get(pg,(float('inf'),))
    extra_cost_tuple = min(old_extra_cost_tuple,extra_cost_tuple)
    pg_to_extra_cost_tuple_dict[pg] = extra_cost_tuple

  #print('MMRQSPXIKX s0_formula_data_list={s0_formula_data_list}'.format(
  #  s0_formula_data_list=s0_formula_data_list
  #))

  for formula_data in s0_formula_data_list:
    print('JHJLGOVYPU formula_data={formula_data}'.format(
      formula_data=formula_data
    ))

  pg_to_min_extra_cost_tuple_dict = {}
  for pg, formula_data_list in pg_to_formula_data_list_dict.items():
    extra_cost_tuple_list = list(map(lambda i: i['extra_cost_tuple'], formula_data_list))
    pg_to_min_extra_cost_tuple_dict[pg] = min(extra_cost_tuple_list)

  s0_formula_data_list = list(filter(
    lambda i: i['extra_cost_tuple']==pg_to_min_extra_cost_tuple_dict[i['product']],
    s0_formula_data_list
  ))

  for formula_data in s0_formula_data_list:
    print('STUZBEMMKR formula_data={formula_data}'.format(
      formula_data=formula_data
    ))

  s1_formula_data_list = []
  needed_gene_set = set()
  for formula_data in s0_formula_data_list:
    pg = formula_data['product']
    if '1' not in pg:
      needed_gene_set.add(pg)
    if pg not in needed_gene_set: continue
    parent_list = formula_data['parent_list']
    if parent_list:
      for g0 in parent_list:
        needed_gene_set.add(g0)
    s1_formula_data_list.append(formula_data)

  for formula_data in s1_formula_data_list:
    print('HVBDCRBNDI formula_data={formula_data}'.format(
      formula_data=formula_data
    ))

  s2_formula_data_list = list(filter(
    lambda i: (i['begin_depth']<=0) and (i['add_depth']>0),
    s1_formula_data_list
  ))

  for formula_data in s2_formula_data_list:
    print('HFSWZFLNRL formula_data={formula_data}'.format(
      formula_data=formula_data
    ))
