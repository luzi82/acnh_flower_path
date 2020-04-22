import futsu.json
import heapq
import math
import os

S6PASS = 0.9999966
S6ERR = 1-S6PASS

s_to_v_list_dict = {
  '00':'0',
  '01':'01',
  '02':'1',
  '11':'0112',
  '12':'12',
  '22':'2',
}

def cross(gene0, gene1):
  v_list_list = []
  for i in range(len(gene0)):
    s = gene0[i] + gene1[i]
    s = ''.join(list(sorted(s)))
    v_list_list.append(s_to_v_list_dict[s])
  
  gene_list = ['']
  for i in range(len(gene0)):
    gene_list_0 = []
    v_list = v_list_list[i]
    for gene in gene_list:
      for v in v_list:
        gene_list_0.append(gene+v)
    gene_list = gene_list_0

  gene_set = set(gene_list)
  
  ret_list = []
  for gene in sorted(gene_set):
    chance = len(list(filter(lambda i:i==gene, gene_list)))
    ret_list.append({'g':gene,'chance':chance})

  chance_sum = len(gene_list)
  for data in ret_list:
    data['p'] = data['chance'] / chance_sum

  return ret_list, chance_sum

def cross_data_list_to_c_to_p_dict(cross_data_list, g_to_c_dict):
  c_to_p_dict = {}
  for cross_data in cross_data_list:
    c = g_to_c_dict[cross_data['g']]
    p = cross_data['p']
    c_to_p_dict[c] = c_to_p_dict.get(c,0) + p
  return c_to_p_dict

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

def verify(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict):
  cross_data_list, chance_sum = cross(parent_gene_0, parent_gene_1)
  cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(cross_data_list, g_to_c_dict)

  product_c_to_cross_data_list_dict = {}
  for cross_c in cross_c_to_p_dict:
    product_c_to_cross_data_list_dict[cross_c] = []
  for cross_data in cross_data_list:
    product_g = cross_data['g']
    product_c = g_to_c_dict[product_g]
    product_c_to_cross_data_list_dict[product_c].append(cross_data)

  for product_c, c_cross_data_list in product_c_to_cross_data_list_dict.items():
    if len(c_cross_data_list) <= 1: continue
    product_gene_list = list(map(lambda i:i['g'], c_cross_data_list))

    verify_g_to_c_set_dict = {}
    verify_gc_to_p_dict = {}
    for product_gene in product_gene_list:
      verify_cross_data_list, _ = cross(verify_gene, product_gene)
      verify_c_to_p_dict = cross_data_list_to_c_to_p_dict(verify_cross_data_list, g_to_c_dict)
      for c, p in verify_c_to_p_dict.items():
        verify_gc_to_p_dict[(product_gene, c)] = p

    for product_gene in product_gene_list:
      # t0
      g_to_p_dict = {
        c_cross_data['g']: c_cross_data['p']/cross_c_to_p_dict[product_c]
        for c_cross_data in c_cross_data_list
      }
      
      # t1-100
      for t in range(100):
        g_to_p_dict0 = {}
        for product_gene0 in product_gene_list:
          p = 0
          TODO

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
  
  for gene_data in gene_data_list:
    if 's' not in gene_data: gene_data['s'] = 0
    if 'm' not in gene_data: gene_data['m'] = 0
  
  g_to_c_dict = {
    i['g']: i['c']
    for i in gene_data_list
  }

  gene_done_set = set()

  depth_gene_heap = []
  for gene_data in gene_data_list:
    if gene_data['s'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))
    if args_myth and gene_data['m'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))

  print(depth_gene_heap)

  depth_gene_list = []
  gene_to_depth_dict = {}
  gene_to_formula_data_list_dict = {}

  def add_formul_data(formula_data):
      heapq.heappush(depth_gene_heap,(formula_data['total_depth'],formula_data['product']))
      if formula_data['product'] not in gene_to_formula_data_list_dict:
        gene_to_formula_data_list_dict[formula_data['product']] = []
      gene_to_formula_data_list_dict[formula_data['product']].append(formula_data)

  for gene_data in gene_data_list:
    if gene_data['s'] > 0:
      add_formul_data({
        'product': gene_data['g'],
        'parent_list': None,
        'step_depth': 0,
        'step_count': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'seed',
      })
    if args_myth and gene_data['m'] > 0:
      add_formul_data({
        'product': gene_data['g'],
        'parent_list': None,
        'step_depth': 0,
        'step_count': 0,
        'add_depth': 0,
        'total_depth': 0,
        'method': 'myth',
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
          'parent_list': (gene, gene_done),
          'step_depth': roll0['step_depth'],
          'step_count': roll0['step_count'],
          'add_depth': add_depth,
          'total_depth': total_depth,
          'method': 'roll',
        })

  print(len(depth_gene_list))
  print(depth_gene_list)
  print(list(sorted(map(lambda i:i[1],depth_gene_list))))

  for depth, gene in depth_gene_list:
    formula_data_list = None
    if gene in gene_to_formula_data_list_dict:
      formula_data_list = gene_to_formula_data_list_dict[gene]
      formula_data_list = filter(lambda i:i['total_depth'] < depth + 0.00001, formula_data_list)
      formula_data_list = list(formula_data_list)
    if formula_data_list:
      for formula_data in formula_data_list:
        if formula_data['total_depth'] > (depth + 0.00001): continue
        print(formula_data)
    else:
      print('{}'.format(gene))
