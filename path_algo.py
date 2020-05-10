from collections import deque
import math

#ACCEPTABLE_ERR = 1/1000000
ACCEPTABLE_ERR = 1/13983816 # mark six
ACCEPTABLE_MISS_P = (1/(1-ACCEPTABLE_ERR))-1

s_to_v_list_dict = {
  '00':'0',
  '01':'01',
  '02':'1',
  '11':'0112',
  '12':'12',
  '22':'2',
}

#def cross_verify(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict, done_depth, product_to_depth_dict):
#  return []

def cross_verify(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict, done_depth, product_to_depth_dict):
  print('FIETHLFSFJ p={parent_gene_0},{parent_gene_1} v={verify_gene}'.format(
    parent_gene_0=parent_gene_0,
    parent_gene_1=parent_gene_1,
    verify_gene=verify_gene,
  ))

  cross_data_list, chance_sum = cross(parent_gene_0, parent_gene_1)
  cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(cross_data_list, g_to_c_dict)

  product_gene_to_cross_data_dict = {
    cross_data['g']: cross_data
    for cross_data in cross_data_list
  }

  #print('UHFUGYMTGG '+str(cross_data_list))

  product_c_to_cross_data_list_dict = {}
  for cross_c in cross_c_to_p_dict:
    product_c_to_cross_data_list_dict[cross_c] = []
  for cross_data in cross_data_list:
    product_g = cross_data['g']
    product_c = g_to_c_dict[product_g]
    product_c_to_cross_data_list_dict[product_c].append(cross_data)

  # return
  cross_verify_data_list = []

  for product_c, c_cross_data_list in product_c_to_cross_data_list_dict.items():
    if len(c_cross_data_list) <= 1: continue
    product_gene_list = list(map(lambda i: i['g'], c_cross_data_list))

    verify_g_to_c_set_dict = {}
    verify_gc_to_logp_dict = {}
    verify_gc_to_p_dict = {}
    verify_g_to_c_to_p_dict_dict = {}
    for product_gene in product_gene_list:
      verify_cross_data_list, _ = cross(verify_gene, product_gene)
      verify_c_to_p_dict = cross_data_list_to_c_to_p_dict(verify_cross_data_list, g_to_c_dict)
      verify_g_to_c_to_p_dict_dict[product_gene] = verify_c_to_p_dict
      for c, p in verify_c_to_p_dict.items():
        verify_gc_to_logp_dict[(product_gene, c)] = math.log2(p)
        verify_gc_to_p_dict[(product_gene, c)] = p
      verify_g_to_c_set_dict[product_gene] = set(verify_c_to_p_dict.keys())

    verify_c_list = set()
    for g, c_set in verify_g_to_c_set_dict.items():
      verify_c_list = verify_c_list | c_set
    verify_c_list = list(sorted(verify_c_list))

    #print('JVJXZYFEMC '+str(verify_g_to_c_set_dict))
    #print('RNVJDDRZRG '+str(verify_gc_to_p_dict))
    #print('MXLUZKAYMJ '+str(verify_gc_to_logp_dict))
    #print('DFOAIHDVDU verify_c_list={verify_c_list}'.format(
    #  verify_c_list=verify_c_list
    #))

    # disable same gene
    diable_gene_set = set()
    for product_gene_0 in product_gene_list:
      for product_gene_1 in product_gene_list:
        if(product_gene_0==product_gene_1):continue
        verify_c_to_p_dict_0 = verify_g_to_c_to_p_dict_dict[product_gene_0]
        verify_c_to_p_dict_1 = verify_g_to_c_to_p_dict_dict[product_gene_1]
        #print('MKFNKVRKLU verify_c_to_p_dict_0={verify_c_to_p_dict_0}'.format(verify_c_to_p_dict_0=verify_c_to_p_dict_0))
        #print('SGYVMCQUTP verify_c_to_p_dict_1={verify_c_to_p_dict_1}'.format(verify_c_to_p_dict_1=verify_c_to_p_dict_1))
        if verify_c_to_p_dict_0==verify_c_to_p_dict_1:
          diable_gene_set.add(product_gene_0)
          diable_gene_set.add(product_gene_1)

    #print('WSHHXUNMFY diable_gene_set={diable_gene_set}'.format(diable_gene_set=diable_gene_set))

    for product_gene in product_gene_list:

      #print('PSWQKIQBMK product_gene='+product_gene)
      if product_gene == parent_gene_0:
        continue
      if product_gene == parent_gene_1:
        continue
      if product_gene == verify_gene:
        continue
      if product_gene in diable_gene_set:
        continue

      verify_gc_to_logp_dict0 = {}
      for c in verify_c_list:
        if (product_gene, c) not in verify_gc_to_logp_dict: continue
        c_logp_base = verify_gc_to_logp_dict.get((product_gene, c),float('-inf'))
        for g0 in product_gene_list:
          verify_gc_to_logp_dict0[(g0,c)] = verify_gc_to_logp_dict.get((g0,c),float('-inf')) - c_logp_base

      #print('HIYQDMQCXS product_gene={} verify_gc_to_logp_dict0={}'.format(product_gene, verify_gc_to_logp_dict0))

      # t0
      g_to_p_dict = {
        c_cross_data['g']: c_cross_data['p']/cross_c_to_p_dict[product_c]
        for c_cross_data in c_cross_data_list
      }
      product_gene_logp_base = math.log2(g_to_p_dict[product_gene])
      g_to_logp_dict = {
        g: math.log2(p) - product_gene_logp_base
        for g, p in g_to_p_dict.items()
      }
      #print('JBGAYVAORB g_to_p_dict='+str(g_to_p_dict))
      #print('CJDPHWVAMT g_to_logp_dict='+str(g_to_logp_dict))
      
      step_limit = product_to_depth_dict.get(product_gene,float('inf')) - done_depth - 1/product_gene_to_cross_data_dict[product_gene]['p']
      step_count = 0
      uncertain_p = 1
      last_c_history_tuple_sum = 0
      step_count_over_break = False
      c_history_tuple_to_p_dict = {}
      c_history_tuple_done = set()
      c_history_tuple_queue = deque()
      
      start_c_history_tuple = tuple([0]*len(verify_c_list))
      c_history_tuple_to_p_dict[start_c_history_tuple] = 1
      c_history_tuple_queue.append(start_c_history_tuple)

      while uncertain_p > ACCEPTABLE_ERR and len(c_history_tuple_queue) > 0:

        c_history_tuple = c_history_tuple_queue.popleft()
        c_history_p = c_history_tuple_to_p_dict[c_history_tuple]

        if c_history_tuple in c_history_tuple_done: continue
        c_history_tuple_done.add(c_history_tuple)

        last_c_history_tuple_sum = sum(c_history_tuple)

        if step_count + uncertain_p * last_c_history_tuple_sum > step_limit:
          step_count_over_break = True
          break

        #print('QDXZRZJJUZ c_history_tuple={c_history_tuple} c_history_p={c_history_p}'.format(
        #  c_history_tuple=c_history_tuple,
        #  c_history_p=c_history_p
        #))

        miss_p = 0
        for product_gene0 in product_gene_list:
          if product_gene0 == product_gene: continue
          logp_list = [
            0 if c_history_tuple[ci] == 0 else
            float('-inf') if (product_gene0,verify_c_list[ci]) not in verify_gc_to_logp_dict else
            c_history_tuple[ci] * verify_gc_to_logp_dict[(product_gene0,verify_c_list[ci])]
            for ci in range(len(verify_c_list))
          ]
          g_p = sum(logp_list)
          g_p = 2**g_p
          g_p *= g_to_p_dict[product_gene0]
          #print('ACZZNHPKJC product_gene0={product_gene0} logp_list={logp_list} g_p={g_p}'.format(
          #  product_gene0=product_gene0,
          #  logp_list=logp_list,
          #  g_p=g_p
          #))
          miss_p += g_p
        #print('XVZFJWVTII miss_p={miss_p}'.format(miss_p=miss_p))

        if miss_p < ACCEPTABLE_MISS_P:
          uncertain_p -= c_history_p
          step_count += c_history_p*sum(c_history_tuple)
          continue

        for ci in range(len(verify_c_list)):
          c = verify_c_list[ci]
          p = verify_gc_to_p_dict.get((product_gene,c),0)
          if p <= 0: continue
          p *= c_history_p
          c_history0_tuple = list(c_history_tuple)
          c_history0_tuple[ci] += 1
          c_history0_tuple = tuple(c_history0_tuple)
          pp = c_history_tuple_to_p_dict.get(c_history0_tuple,0)
          c_history_tuple_to_p_dict[c_history0_tuple] = pp+p
          c_history_tuple_queue.append(c_history0_tuple)

      if step_count_over_break: continue
      if uncertain_p > ACCEPTABLE_ERR: continue

      add_step = step_count + 1/product_gene_to_cross_data_dict[product_gene]['p']
      add_step += uncertain_p * last_c_history_tuple_sum
      cross_verify_data_list.append({
        'product': product_gene,
        'product.color': g_to_c_dict[product_gene],
        'parent_list': (parent_gene_0,parent_gene_1,verify_gene),
        'step_depth': 1,
        'step_count': add_step,
        'add_depth': add_step,
        'method': 'cross_verify',
      })

  return cross_verify_data_list

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
