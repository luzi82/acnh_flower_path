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
  gene_to_formula_list_dict = {}

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

      #print('---')
      #print(gene)
      #print(gene_done)
      #print(c_to_gene_set_dict)
      #print('---')

      for gene_chance in gene_chance_list:
        g = gene_chance['g']
        if g in gene_done_set: continue
        c = g_to_c_dict[g]
        gene_set = c_to_gene_set_dict[c]
        if len(gene_set) > 1: continue
        chance = gene_chance['chance']
        #depth0 = math.log(chance_sum/chance, 2)
        depth0 = chance_sum/chance
        depth1 = depth0 + depth
        
        heapq.heappush(depth_gene_heap,(depth1,g))

        if g not in gene_to_formula_list_dict:
          gene_to_formula_list_dict[g] = []
        
        gene_to_formula_list_dict[g].append((gene_done, gene, depth0, depth1))

        #print('---')
        #print(gene     +' '+g_to_c_dict[gene])
        #print(gene_done+' '+g_to_c_dict[gene_done])
        #print(g        +' '+g_to_c_dict[g])
        #print('---')

  print(depth_gene_list)
  #print(gene_to_formula_list_dict)

  for depth, gene in depth_gene_list:
    formula_list = None
    if gene in gene_to_formula_list_dict:
      formula_list = gene_to_formula_list_dict[gene]
      formula_list = filter(lambda i:i[3] < depth + 0.00001, formula_list)
      formula_list = list(formula_list)
    if formula_list:
      for gene0, gene1, depth0, depth1 in formula_list:
        if depth1 > depth + 0.00001: continue
        print('{} {} {} {} {}'.format(gene0, gene1, gene, depth0, depth1))
    else:
      print('{}'.format(gene))
