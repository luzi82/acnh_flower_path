import run
import futsu.json
import os

def ssorted(ddict):
  return {k:ddict[k] for k in sorted(ddict.keys())}

if __name__ == '__main__':

  verify_data_path = 'verify.json'
  verify_data = futsu.json.path_to_data(verify_data_path)

  flower_data_path = os.path.join('res','{}.json'.format(verify_data['flower_name']))
  flower_data = futsu.json.path_to_data(flower_data_path)

  gene_data_list = flower_data['gene_data_list']
  
  for gene_data in gene_data_list:
    if 's' not in gene_data: gene_data['s'] = 0
    if 'm' not in gene_data: gene_data['m'] = 0
  
  g_to_c_dict = {
    i['g']: i['c']
    for i in gene_data_list
  }

  pool_data_list = verify_data['pool_data_list']
  pool_gene_list = list(map(lambda i:i['g'], pool_data_list))

  target_gene_list = filter(lambda i:i.get('t',0)>0, pool_data_list)
  target_gene_list = list(map(lambda i:i['g'], target_gene_list))

  own_gene_list = verify_data['own_gene_list']

  for own_gene in own_gene_list:
    pool_g_to_c_to_p_dict_dict = {}
    for pool_gene in pool_gene_list:
      cross_data_list, _ = run.cross(own_gene, pool_gene)
      c_to_p_dict = run.cross_data_list_to_c_to_p_dict(cross_data_list, g_to_c_dict)
      pool_g_to_c_to_p_dict_dict[pool_gene] = ssorted(c_to_p_dict)

    print(own_gene)
    print(pool_g_to_c_to_p_dict_dict)

    target_c_remain_set = set()
    for target_gene in target_gene_list:
      tmp = set(pool_g_to_c_to_p_dict_dict[target_gene].keys())
      target_c_remain_set = target_c_remain_set | tmp

    #print(target_c_remain_set)
    for pool_gene, c_to_p_dict in pool_g_to_c_to_p_dict_dict.items():
      if pool_gene in target_gene_list: continue
      c_set = set(c_to_p_dict.keys())
      #print(c_set)
      target_c_remain_set = target_c_remain_set - c_set
      #print(target_c_remain_set)
    
    if len(target_c_remain_set) <= 0:
      continue
    
    for target_gene in target_gene_list:
      p_total = pool_g_to_c_to_p_dict_dict[target_gene]
      p_total = filter(lambda i: i[0] in target_c_remain_set , p_total.items())
      p_total = list(p_total)
      
      color_list = list(map(lambda i: i[0], p_total))
      
      p_total = map(lambda i: i[1], p_total)
      p_total = sum(p_total)
      
      print('{target_gene} {own_gene} {color_list} {p_total}'.format(
        target_gene = target_gene,
        own_gene = own_gene,
        color_list = str(color_list),
        p_total = p_total
      ))
