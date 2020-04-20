import futsu.json
import heapq
import math

F='res/tulip.json'

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

  return ret_list, len(gene_list)

if __name__ == '__main__':

  flower_data = futsu.json.path_to_data(F)
  
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
    if gene_data['m'] > 0:
      heapq.heappush(depth_gene_heap,(0,gene_data['g']))

  print(depth_gene_heap)

  depth_gene_list = []
  gene_to_depth_dict = {}
  gene_to_formula_data_list_dict = {}

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
        #depth0 = math.log(chance_sum/chance, 2)
        add_depth = chance_sum/chance
        total_depth = add_depth + depth
        
        heapq.heappush(depth_gene_heap,(total_depth,g))

        if g not in gene_to_formula_data_list_dict:
          gene_to_formula_data_list_dict[g] = []
        
        gene_to_formula_data_list_dict[g].append({
          'parent_list': (gene_done, gene),
          'add_depth': add_depth,
          'total_depth': total_depth,
        })

  print(depth_gene_list)

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
