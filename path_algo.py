import heapq
import math
import numpy.linalg
import numpy as np

from collections import deque
from fractions import Fraction


#ACCEPTABLE_ERR = 1/1000000
ACCEPTABLE_ERR = 1/13983816 # mark six
ACCEPTABLE_MISS_P = (1/(1-ACCEPTABLE_ERR))-1
FLOAT_CORRECT = 0.00001

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
  #print('FIETHLFSFJ p={parent_gene_0},{parent_gene_1} v={verify_gene}'.format(
  #  parent_gene_0=parent_gene_0,
  #  parent_gene_1=parent_gene_1,
  #  verify_gene=verify_gene,
  #))

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

def roll(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict, done_depth, product_to_depth_dict):
  print('LPSDLVHKMM p={parent_gene_0},{parent_gene_1} v={verify_gene}'.format(
    parent_gene_0=parent_gene_0,
    parent_gene_1=parent_gene_1,
    verify_gene=verify_gene,
  ))

  cross0_data_list, cross0_chance_sum = cross(parent_gene_0, parent_gene_1)
  cross0_c_to_p_dict = cross_data_list_to_c_to_p_dict(cross0_data_list, g_to_c_dict)

  g0_to_cross_data_dict = {
    cross_data['g']: cross_data
    for cross_data in cross0_data_list
  }

  c0_to_cross_data_list_dict = {}
  for c0 in cross0_c_to_p_dict:
    c0_to_cross_data_list_dict[c0] = []
  for cross_data in cross0_data_list:
    g0 = cross_data['g']
    c0 = g_to_c_dict[g0]
    c0_to_cross_data_list_dict[c0].append(cross_data)

  roll_data_list = []

  for c0, c0_cross_data_list in c0_to_cross_data_list_dict.items():
    if len(c0_cross_data_list) != 2: continue
    print('DGYWTCWEXD c0={}, c0_cross_data_list={}'.format(c0,c0_cross_data_list))

    g0_list = list(map(lambda i: i['g'], c0_cross_data_list))
    g0_set = set(g0_list)

    # check merge
    exist_one = False
    exist_all = False
    bad = False
    tg = None # target gene
    for g0 in g0_list:
      cross1_data_list, _ = cross(verify_gene, g0)
      g1_itr = map(lambda i:i['g'], cross1_data_list)
      c0g1_set = set(filter(lambda i:g_to_c_dict[i]==c0, g1_itr))
      print('DKYWEPEEJO g0_set={g0_set}, g0={g0} c0g1_set={c0g1_set}'.format(
        g0_set=g0_set,
        g0=g0,
        c0g1_set=c0g1_set,
      ))
      if len(c0g1_set) == 2:
        if c0g1_set != g0_set:
          bad = True
          break
        exist_all = True
      elif len(c0g1_set) == 1:
        tg = list(c0g1_set)[0]
        if tg not in g0_set:
          bad = True
          break
        exist_one = True
      else:
        bad = True
        break
    if bad: continue
    if not exist_one: continue
    if not exist_all: continue
    
    fg = list(g0_set-set([tg]))[0] # from gene

    print('AWXXOJNZMV fg={}, tg={}'.format(fg,tg))

    step_limit = product_to_depth_dict.get(tg,float('inf')) - done_depth - 1/cross0_c_to_p_dict[c0]
    if step_limit < 1: continue

    fg_cross_data_list, _ = cross(verify_gene, fg)
    fg_cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(fg_cross_data_list, g_to_c_dict)
    fg_cross_c_to_pf_dict = {
      c: Fraction(p)
      for c, p in fg_cross_c_to_p_dict.items()
    }
    fg_fg_pf = Fraction(list(filter(lambda i:i['g']==fg,fg_cross_data_list))[0]['p'])
    fg_tg_pf = Fraction(list(filter(lambda i:i['g']==tg,fg_cross_data_list))[0]['p'])

    tg_cross_data_list, _ = cross(verify_gene, tg)
    tg_cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(tg_cross_data_list, g_to_c_dict)
    tg_cross_c_to_pf_dict = {
      c: Fraction(p)
      for c, p in tg_cross_c_to_p_dict.items()
    }

    fg_state_ffpredf_ftpref_ttpref_p_list = [
      ((False, fg_cross_c_to_pf_dict.get(c,0), Fraction(0), tg_cross_c_to_pf_dict.get(c,0)), float(fg_cross_c_to_pf_dict[c]))
      for c in fg_cross_c_to_pf_dict if c != c0
    ]
    fg_state_ffpredf_ftpref_ttpref_p_list.append(
      ((False, fg_fg_pf, fg_tg_pf, tg_cross_c_to_pf_dict.get(c0,0)), float(fg_fg_pf))
    )
    fg_state_ffpredf_ftpref_ttpref_p_list.append(
      ((True,  fg_fg_pf, fg_tg_pf, tg_cross_c_to_pf_dict.get(c0,0)), float(fg_tg_pf))
    )
    fg_state_ffpredf_ftpref_ttpref_to_p_dict = {}
    for state_ffpredf_ftpref_ttpref, p in fg_state_ffpredf_ftpref_ttpref_p_list:
      pp = fg_state_ffpredf_ftpref_ttpref_to_p_dict.get(state_ffpredf_ftpref_ttpref, 0)
      fg_state_ffpredf_ftpref_ttpref_to_p_dict[state_ffpredf_ftpref_ttpref] = pp + p
    fg_state_ffpredf_ftpref_ttpref_p_list = [
      (state_ffpredf_ftpref_ttpref, p)
      for state_ffpredf_ftpref_ttpref, p in fg_state_ffpredf_ftpref_ttpref_to_p_dict.items()
    ]

    tg_state_ffpredf_ftpref_ttpref_p_list = [
      ((True, fg_cross_c_to_pf_dict.get(c,0), Fraction(0), tg_cross_c_to_pf_dict.get(c,0)), float(tg_cross_c_to_pf_dict[c]))
      for c in tg_cross_c_to_pf_dict if c != c0
    ]
    tg_state_ffpredf_ftpref_ttpref_p_list.append(
      ((True,  fg_fg_pf, fg_tg_pf, tg_cross_c_to_pf_dict.get(c0,0)), float(tg_cross_c_to_pf_dict[c0]))
    )
    tg_state_ffpredf_ftpref_ttpref_to_p_dict = {}
    for state_ffpredf_ftpref_ttpref, p in tg_state_ffpredf_ftpref_ttpref_p_list:
      pp = tg_state_ffpredf_ftpref_ttpref_to_p_dict.get(state_ffpredf_ftpref_ttpref, 0)
      tg_state_ffpredf_ftpref_ttpref_to_p_dict[state_ffpredf_ftpref_ttpref] = pp + p
    tg_state_ffpredf_ftpref_ttpref_p_list = [
      (state_ffpredf_ftpref_ttpref, p)
      for state_ffpredf_ftpref_ttpref, p in tg_state_ffpredf_ftpref_ttpref_to_p_dict.items()
    ]
    
    state_to_state_ffpredf_ftpref_ttpref_p_list_dict = {
      False: fg_state_ffpredf_ftpref_ttpref_p_list,
      True:  tg_state_ffpredf_ftpref_ttpref_p_list,
    }

    uncertain_p = 1

    p_state_predf_heap = []
    state_predf_to_p_dict = {}
    state_predf_to_mat_idx_dict = {}
    mat_idx_to_state_predf_dict = {}
    state_predf_to_formula_data_dict = {}

    def add_state_predf_p(state_predf, p):
      nonlocal uncertain_p
      if state_predf in state_predf_to_formula_data_dict:
        uncertain_p -= p
      else:
        pp = p + state_predf_to_p_dict.get(state_predf, 0)
        state_predf_to_p_dict[state_predf] = pp
        heapq.heappush(p_state_predf_heap, (-pp, state_predf))
      if state_predf not in state_predf_to_mat_idx_dict:
        mat_idx = len(state_predf_to_mat_idx_dict)
        state_predf_to_mat_idx_dict[state_predf] = mat_idx
        mat_idx_to_state_predf_dict[mat_idx] = state_predf
    
    predf = Fraction(
      g0_to_cross_data_dict[tg]['chance'],
      g0_to_cross_data_dict[tg]['chance']+g0_to_cross_data_dict[fg]['chance']
    )

    i0_state_predf = (None, None)
    i0_formula_data = {
      'state_predf_p_list':[
        ((True, predf), float(predf)),
        ((False, predf), 1-float(predf)),
      ],
      'else': 0
    }
    state_predf_to_mat_idx_dict[i0_state_predf] = 0
    mat_idx_to_state_predf_dict[0] = i0_state_predf
    state_predf_to_formula_data_dict[i0_state_predf] = i0_formula_data

    add_state_predf_p((True, predf), float(predf))
    add_state_predf_p((False, predf), 1-float(predf))

    while((len(p_state_predf_heap)>0)and(uncertain_p>ACCEPTABLE_ERR)):
      _, state_predf = heapq.heappop(p_state_predf_heap)
      if state_predf in state_predf_to_formula_data_dict:
        continue
      state, predf = state_predf
      state_predf__p = state_predf_to_p_dict[state_predf]
      print('XQKNNLQJTN state_predf={state_predf} state_predf__p={state_predf__p} uncertain_p={uncertain_p}'.format(
        state_predf=state_predf,
        state_predf__p=state_predf__p,
        uncertain_p=uncertain_p,
      ))
      if state and (float(predf) > 1- ACCEPTABLE_ERR):
        state_predf_to_formula_data_dict[state_predf] = {
          'state_predf_p_list':[],
          'else': 0
        }
        uncertain_p -= state_predf__p
        continue

      new_state_predf_p_list = []
      state_ffpredf_ftpref_ttpref_p_list = state_to_state_ffpredf_ftpref_ttpref_p_list_dict[state]
      for state_ffpredf_ftpref_ttpref, new_p in state_ffpredf_ftpref_ttpref_p_list:
        new_state, new_ffpredf, new_ftpref, new_ttpref = state_ffpredf_ftpref_ttpref
        new_fpredf = (1-predf) * new_ffpredf
        new_tpredf = ((1-predf) * new_ftpref) + (predf * new_ttpref)
        new_predf  = new_tpredf / (new_fpredf+new_tpredf)
        new_state_predf_p_list.append(
          ((new_state, new_predf), new_p)
        )

      state_predf_to_formula_data_dict[state_predf] = {
        'state_predf_p_list':new_state_predf_p_list,
        'else': 1
      }

      print('RMBHBLCAAI new_state_predf_p_list={new_state_predf_p_list}'.format(
        new_state_predf_p_list=new_state_predf_p_list
      ))

      for new_state_predf, new_p in new_state_predf_p_list:
        add_state_predf_p(new_state_predf, new_p*state_predf__p)

      # cal over
      grid_size = len(state_predf_to_mat_idx_dict)
      a_ll_np = np.zeros((grid_size,grid_size))
      b_l_np = np.zeros((grid_size,))
      for idx0 in range(len(state_predf_to_mat_idx_dict)):
        a_ll_np[idx0,idx0] = 1
        state_predf_0 = mat_idx_to_state_predf_dict[idx0]
        if state_predf_0 not in state_predf_to_formula_data_dict: continue
        formula_data_0 = state_predf_to_formula_data_dict[state_predf_0]
        for state_predf_1, p1 in formula_data_0['state_predf_p_list']:
          idx1 = state_predf_to_mat_idx_dict[state_predf_1]
          a_ll_np[idx0,idx1] -= p1
        b_l_np[idx0] = formula_data_0['else']
      x_l_np = np.linalg.solve(a_ll_np, b_l_np)
      add_depth = x_l_np[0]
      if add_depth - step_limit > FLOAT_CORRECT:
        bad = True
        break

    if bad: continue
    
    grid_size = len(state_predf_to_formula_data_dict)
    a_ll_np = np.zeros((grid_size,grid_size))
    b_l_np = np.zeros((grid_size,))
    for idx0 in range(len(state_predf_to_mat_idx_dict)):
      a_ll_np[idx0,idx0] = 1
      state_predf_0 = mat_idx_to_state_predf_dict[idx0]
      if state_predf_0 not in state_predf_to_formula_data_dict:
        a_ll_np[idx0,0] -= 1
      else:
        formula_data_0 = state_predf_to_formula_data_dict[state_predf_0]
        for state_predf_1, p1 in formula_data_0['state_predf_p_list']:
          idx1 = state_predf_to_mat_idx_dict[state_predf_1]
          a_ll_np[idx0,idx1] -= p1
        b_l_np[idx0] = formula_data_0['else']
    x_l_np = np.linalg.solve(a_ll_np, b_l_np)
    add_depth = x_l_np[0]
    if add_depth - step_limit > FLOAT_CORRECT:
      break

    roll_data_list.append({
      'product': tg,
      'product.color': g_to_c_dict[tg],
      'parent_list': (parent_gene_0,parent_gene_1,verify_gene),
      'step_depth': 1,
      'step_count': add_step,
      'add_depth': add_step,
      'method': 'roll',
    })

  print('SPEMSCKDFN roll_data_list='+str(roll_data_list))

  return roll_data_list

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
