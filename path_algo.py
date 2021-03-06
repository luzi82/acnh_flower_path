import copy
import collections
import heapq
import json
import math
import numpy.linalg
import numpy as np

from collections import deque
#from fractions import Fraction

def ll(v):
  if v == 0: return float('-inf')
  return math.log2(v)*LOG2_PRECISION

def rr(v):
  if v == float('-inf'): return v
  return round(v)

def kk(v):
  return 2**(v/8)

#ACCEPTABLE_ERR = 1/1000000
ACCEPTABLE_ERR = 1/13983816 # mark six
LOG2_PRECISION = 8
ACCEPTABLE_ERR_LOGS2 = ll(ACCEPTABLE_ERR)
ACCEPTABLE_MISS_P = (1/(1-ACCEPTABLE_ERR))-1
FLOAT_CORRECT = 0.00001
BIG_MISS_GIVEIP_MISS_P = 0.9
BIG_MISS_STEP_COUNT = 20
BIG_MISS_GIVEIP_MISS_P0 = 0.5

EARLY_GIVEUP_STEP_COUNT = 10
EARLY_GIVEUP_BAD_P = 0.75

#print('WVFDOGKNWE ACCEPTABLE_ERR_LOGS2={ACCEPTABLE_ERR_LOGS2}'.format(
#  ACCEPTABLE_ERR_LOGS2=ACCEPTABLE_ERR_LOGS2
#))

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

def cross_verify(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict, done_depth, product_to_depth_dict, path_algo_cache):
  #print('FIETHLFSFJ p={parent_gene_0},{parent_gene_1} v={verify_gene}'.format(
  #  parent_gene_0=parent_gene_0,
  #  parent_gene_1=parent_gene_1,
  #  verify_gene=verify_gene,
  #))

  cross_data_list, chance_sum = cross(parent_gene_0, parent_gene_1, path_algo_cache)
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
      verify_cross_data_list, _ = cross(verify_gene, product_gene, path_algo_cache)
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

    # disable same c_to_p dict gene
    disable_gene_set = set()
    for product_gene_0 in product_gene_list:
      for product_gene_1 in product_gene_list:
        if(product_gene_0==product_gene_1):continue
        if(g_to_c_dict[product_gene_0]!=g_to_c_dict[product_gene_1]):continue
        verify_c_to_p_dict_0 = verify_g_to_c_to_p_dict_dict[product_gene_0]
        verify_c_to_p_dict_1 = verify_g_to_c_to_p_dict_dict[product_gene_1]
        #print('MKFNKVRKLU verify_c_to_p_dict_0={verify_c_to_p_dict_0}'.format(verify_c_to_p_dict_0=verify_c_to_p_dict_0))
        #print('SGYVMCQUTP verify_c_to_p_dict_1={verify_c_to_p_dict_1}'.format(verify_c_to_p_dict_1=verify_c_to_p_dict_1))
        #if verify_c_to_p_dict_0==verify_c_to_p_dict_1:
        #  disable_gene_set.add(product_gene_0)
        #  disable_gene_set.add(product_gene_1)
        verify_c_set_0 = set(verify_c_to_p_dict_0.keys())
        verify_c_set_1 = set(verify_c_to_p_dict_1.keys())
        if verify_c_set_0==verify_c_set_1:
          disable_gene_set.add(product_gene_0)
          disable_gene_set.add(product_gene_1)

    #print('WSHHXUNMFY disable_gene_set={disable_gene_set}'.format(disable_gene_set=disable_gene_set))

    for product_gene in product_gene_list:

      #print('PSWQKIQBMK product_gene='+product_gene)

      if product_gene == parent_gene_0:
        continue
      if product_gene == parent_gene_1:
        continue
      if product_gene == verify_gene:
        continue
      if product_gene in disable_gene_set:
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
      
      #step_limit = product_to_depth_dict.get(product_gene,float('inf')) - done_depth - 1/product_gene_to_cross_data_dict[product_gene]['p']
      add_step_limit = product_to_depth_dict.get(product_gene,float('inf')) - done_depth
      step_count = 0
      uncertain_p = 1
      big_hit_p = 0
      big_miss_p = 0
      big_miss_ps = 0
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

        if big_miss_p > BIG_MISS_GIVEIP_MISS_P0:
          #print('FSMREWMMLA break')
          break

        if last_c_history_tuple_sum >= BIG_MISS_STEP_COUNT:
          #print('KLMMIWJGHZ break')
          big_miss_p += uncertain_p
          big_miss_ps += uncertain_p*last_c_history_tuple_sum
          uncertain_p = 0
          break

        if last_c_history_tuple_sum >= EARLY_GIVEUP_STEP_COUNT and uncertain_p+big_miss_p >= EARLY_GIVEUP_BAD_P:
          #print('MRPRLWGROP break')
          break

        tmp_add_step = 0 if big_hit_p == 0 else step_count / big_hit_p
        tmp_add_step = tmp_add_step + 1/product_gene_to_cross_data_dict[product_gene]['p']
        tmp_add_step += last_c_history_tuple_sum * uncertain_p
        tmp_add_step = ((1-big_miss_p)*tmp_add_step + (big_miss_ps))/(1-big_miss_p)
        
        #print('BXSPMQCWBJ step_count={step_count} uncertain_p={uncertain_p} last_c_history_tuple_sum={last_c_history_tuple_sum} tmp_add_step={tmp_add_step} big_miss_p={big_miss_p}'.format(
        #  step_count=step_count,
        #  uncertain_p=uncertain_p,
        #  last_c_history_tuple_sum=last_c_history_tuple_sum,
        #  tmp_add_step=tmp_add_step,
        #  big_miss_p=big_miss_p
        #))

        if tmp_add_step - add_step_limit > FLOAT_CORRECT:
          step_count_over_break = True
          break

        #print('QDXZRZJJUZ c_history_tuple={c_history_tuple} c_history_p={c_history_p} uncertain_p={uncertain_p} big_miss_p={big_miss_p}'.format(
        #  uncertain_p=uncertain_p,
        #  c_history_tuple=c_history_tuple,
        #  c_history_p=c_history_p,
        #  big_miss_p=big_miss_p
        #))

        hit_p = 0
        miss_p = 0
        for product_gene0 in product_gene_list:
          #for ci in range(len(verify_c_list)):
          #  print('HEBIUJFFPQ ci={ci} cht={cht} vgtld={vgtld}'.format(
          #    ci=ci,
          #    cht=c_history_tuple[ci],
          #    vgtld=verify_gc_to_logp_dict0.get((product_gene0,verify_c_list[ci]),float('-inf'))
          #  ))
          logp_list = [
            0 if c_history_tuple[ci] == 0 else
            float('-inf') if (product_gene0,verify_c_list[ci]) not in verify_gc_to_logp_dict0 else
            (c_history_tuple[ci] * verify_gc_to_logp_dict0[(product_gene0,verify_c_list[ci])])
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
          if product_gene0 == product_gene:
            hit_p += g_p
          else:
            miss_p += g_p
        miss_p = miss_p/(miss_p+hit_p)
        #print('XVZFJWVTII miss_p={miss_p}'.format(miss_p=miss_p))

        if miss_p < ACCEPTABLE_MISS_P:
          #print('PBRHHHSUDR miss_p={miss_p}'.format(miss_p=miss_p))
          uncertain_p -= c_history_p
          big_hit_p += c_history_p
          step_count += c_history_p*sum(c_history_tuple)
          continue

        if miss_p > BIG_MISS_GIVEIP_MISS_P:
          uncertain_p -= c_history_p
          big_miss_p  += c_history_p
          big_miss_ps += c_history_p*sum(c_history_tuple)
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

      if step_count_over_break:
        #print('JVQSOWATQU step_count_over_break product_gene={product_gene}'.format(
        #  product_gene=product_gene
        #))
        continue
      if uncertain_p > ACCEPTABLE_ERR:
        #print('RKOTERQIPG uncertain_p product_gene={product_gene}'.format(
        #  product_gene=product_gene
        #))
        continue
      if big_miss_p > BIG_MISS_GIVEIP_MISS_P0:
        #print('GXWPDFNVHI big_miss_p product_gene={product_gene}'.format(
        #  product_gene=product_gene
        #))
        continue

      tmp_add_step = 0 if big_hit_p == 0 else step_count / big_hit_p
      tmp_add_step = tmp_add_step + 1/product_gene_to_cross_data_dict[product_gene]['p']
      tmp_add_step += last_c_history_tuple_sum * uncertain_p
      tmp_add_step = ((1-big_miss_p)*tmp_add_step + (big_miss_ps))/(1-big_miss_p)

      #print('HFPUQARLDG {}'.format(pj({
      #  'step_count' :step_count,
      #  'last_c_history_tuple_sum': last_c_history_tuple_sum,
      #  'uncertain_p':uncertain_p,
      #  'big_hit_p': big_hit_p,
      #  'big_miss_p': big_miss_p,
      #})))

      add_step = tmp_add_step
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

def roll(parent_gene_0, parent_gene_1, verify_gene, g_to_c_dict, done_depth, product_to_depth_dict, path_algo_cache):
  #print('LPSDLVHKMM p={parent_gene_0},{parent_gene_1} v={verify_gene}'.format(
  #  parent_gene_0=parent_gene_0,
  #  parent_gene_1=parent_gene_1,
  #  verify_gene=verify_gene,
  #))

  cross0_data_list, cross0_chance_sum = cross(parent_gene_0, parent_gene_1, path_algo_cache)
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
    #print('DGYWTCWEXD c0={}, c0_cross_data_list={}'.format(c0,c0_cross_data_list))

    g0_list = list(map(lambda i: i['g'], c0_cross_data_list))
    g0_set = set(g0_list)

    # check merge
    exist_one = False
    exist_all = False
    bad = False
    tg = None # target gene
    for g0 in g0_list:
      cross1_data_list, _ = cross(verify_gene, g0, path_algo_cache)
      g1_itr = map(lambda i:i['g'], cross1_data_list)
      c0g1_set = set(filter(lambda i:g_to_c_dict[i]==c0, g1_itr))
      #print('DKYWEPEEJO g0_set={g0_set}, g0={g0} c0g1_set={c0g1_set}'.format(
      #  g0_set=g0_set,
      #  g0=g0,
      #  c0g1_set=c0g1_set,
      #))
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

    #print('AWXXOJNZMV fg={}, tg={}'.format(fg,tg))

    step_limit = product_to_depth_dict.get(tg,float('inf')) - done_depth - 1/cross0_c_to_p_dict[c0]
    if step_limit < 1: continue

    fg_cross_data_list, _ = cross(verify_gene, fg, path_algo_cache)
    fg_cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(fg_cross_data_list, g_to_c_dict)
    fg_cross_c_to_plog_dict = {k: ll(v) for k,v in fg_cross_c_to_p_dict.items()}
    fg_fg_p = list(filter(lambda i:i['g']==fg,fg_cross_data_list))[0]['p']
    fg_fg_plog = ll(fg_fg_p)
    fg_tg_p = list(filter(lambda i:i['g']==tg,fg_cross_data_list))[0]['p']
    fg_tg_plog = ll(fg_tg_p)

    tg_cross_data_list, _ = cross(verify_gene, tg, path_algo_cache)
    tg_cross_c_to_p_dict = cross_data_list_to_c_to_p_dict(tg_cross_data_list, g_to_c_dict)
    tg_cross_c_to_plog_dict = {k: ll(v) for k,v in tg_cross_c_to_p_dict.items()}

    fg_state_ffpredlog_ftpredlog_ttpredlog_p_list = [
      ((False, fg_cross_c_to_plog_dict.get(c,float('-inf')), float('-inf'), tg_cross_c_to_plog_dict.get(c,float('-inf'))), fg_cross_c_to_p_dict[c])
      for c in fg_cross_c_to_plog_dict if c != c0
    ]
    fg_state_ffpredlog_ftpredlog_ttpredlog_p_list.append(
      ((False, fg_fg_plog, fg_tg_plog, tg_cross_c_to_plog_dict.get(c0,float('-inf'))), fg_fg_p)
    )
    fg_state_ffpredlog_ftpredlog_ttpredlog_p_list.append(
      ((True,  fg_fg_plog, fg_tg_plog, tg_cross_c_to_plog_dict.get(c0,float('-inf'))), fg_tg_p)
    )
    fg_state_ffpredlog_ftpredlog_ttpredlog_p_dict = {}
    for state_ffpredlog_ftpredlog_ttpredlog, p in fg_state_ffpredlog_ftpredlog_ttpredlog_p_list:
      pp = fg_state_ffpredlog_ftpredlog_ttpredlog_p_dict.get(state_ffpredlog_ftpredlog_ttpredlog, 0)
      fg_state_ffpredlog_ftpredlog_ttpredlog_p_dict[state_ffpredlog_ftpredlog_ttpredlog] = pp + p
    fg_state_ffpredlog_ftpredlog_ttpredlog_p_list = [
      (state_ffpredlog_ftpredlog_ttpredlog, p)
      for state_ffpredlog_ftpredlog_ttpredlog, p in fg_state_ffpredlog_ftpredlog_ttpredlog_p_dict.items()
    ]

    tg_state_ffpredlog_ftpredlog_ttpredlog_p_list = [
      ((True, fg_cross_c_to_plog_dict.get(c,float('-inf')), float('-inf'), tg_cross_c_to_plog_dict.get(c,float('-inf'))), tg_cross_c_to_p_dict[c])
      for c in tg_cross_c_to_plog_dict if c != c0
    ]
    tg_state_ffpredlog_ftpredlog_ttpredlog_p_list.append(
      ((True,  fg_fg_plog, fg_tg_plog, tg_cross_c_to_plog_dict.get(c0,float('-inf'))), tg_cross_c_to_p_dict[c0])
    )
    tg_state_ffpredlog_ftpredlog_ttpredlog_p_dict = {}
    for state_ffpredlog_ftpredlog_ttpredlog, p in tg_state_ffpredlog_ftpredlog_ttpredlog_p_list:
      pp = tg_state_ffpredlog_ftpredlog_ttpredlog_p_dict.get(state_ffpredlog_ftpredlog_ttpredlog, 0)
      tg_state_ffpredlog_ftpredlog_ttpredlog_p_dict[state_ffpredlog_ftpredlog_ttpredlog] = pp + p
    tg_state_ffpredlog_ftpredlog_ttpredlog_p_list = [
      (state_ffpredlog_ftpredlog_ttpredlog, p)
      for state_ffpredlog_ftpredlog_ttpredlog, p in tg_state_ffpredlog_ftpredlog_ttpredlog_p_dict.items()
    ]
    
    state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict = {
      False: fg_state_ffpredlog_ftpredlog_ttpredlog_p_list,
      True:  tg_state_ffpredlog_ftpredlog_ttpredlog_p_list,
    }
    #print('YUKGBLJRIA state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict={state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict}'.format(
    #  state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict=state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict
    #))

    uncertain_p = 1

    p_state_predlogr_heap = []
    state_predlogr_to_p_dict = {}
    state_predlogr_to_mat_idx_dict = {}
    mat_idx_to_state_predlogr_dict = {}
    state_predlogr_to_formula_data_dict = {}

    def add_state_predlogr_p(state_predlogr, p):
      nonlocal uncertain_p
      if state_predlogr in state_predlogr_to_formula_data_dict:
        uncertain_p -= p
      else:
        pp = p + state_predlogr_to_p_dict.get(state_predlogr, 0)
        state_predlogr_to_p_dict[state_predlogr] = pp
        heapq.heappush(p_state_predlogr_heap, (-pp, state_predlogr))
      if state_predlogr not in state_predlogr_to_mat_idx_dict:
        mat_idx = len(state_predlogr_to_mat_idx_dict)
        state_predlogr_to_mat_idx_dict[state_predlogr] = mat_idx
        mat_idx_to_state_predlogr_dict[mat_idx] = state_predlogr
        #print('FWGGFIUDAD state_predlogr={state_predlogr} mat_idx={mat_idx}'.format(
        #  state_predlogr=state_predlogr,
        #  mat_idx=mat_idx,
        #))
    
    pred = g0_to_cross_data_dict[fg]['chance']/(g0_to_cross_data_dict[tg]['chance']+g0_to_cross_data_dict[fg]['chance'])

    i0_state_predlogr = (None, None)
    i0_formula_data = {
      'state_predlogr_p_list':[
        ((True,  rr(ll(pred))), 1-pred),
        ((False, rr(ll(pred))),   pred),
      ],
      'else': 0
    }
    state_predlogr_to_mat_idx_dict[i0_state_predlogr] = 0
    mat_idx_to_state_predlogr_dict[0] = i0_state_predlogr
    state_predlogr_to_formula_data_dict[i0_state_predlogr] = i0_formula_data

    for i in i0_formula_data['state_predlogr_p_list']:
      add_state_predlogr_p(i[0], i[1])

    while((len(p_state_predlogr_heap)>0)and(uncertain_p>ACCEPTABLE_ERR)):
      _, state_predlogr = heapq.heappop(p_state_predlogr_heap)
      if state_predlogr in state_predlogr_to_formula_data_dict:
        continue
      state, predlogr = state_predlogr
      state_predlogr__p = state_predlogr_to_p_dict[state_predlogr]
      #print('XQKNNLQJTN state_predlogr={state_predlogr} state_predlogr__p={state_predlogr__p} uncertain_p={uncertain_p}'.format(
      #  state_predlogr=state_predlogr,
      #  state_predlogr__p=state_predlogr__p,
      #  uncertain_p=uncertain_p,
      #))

      if predlogr < ACCEPTABLE_ERR_LOGS2:
        if state:
          state_predlogr_to_formula_data_dict[state_predlogr] = {
            'state_predlogr_p_list':[],
            'else': 0
          }
          uncertain_p -= state_predlogr__p
        else:
          state_predlogr_to_formula_data_dict[state_predlogr] = {
            'state_predlogr_p_list':[
              ((False, 0), 1)
            ],
            'else': 0
          }
          add_state_predlogr_p((False, 0),state_predlogr__p)
        continue

      new_state_predlogr_p_list = []
      state_ffpredlog_ftpredlog_ttpredlog_p_list = state_to_state_ffpredlog_ftpredlog_ttpredlog_p_list_dict[state]
      for state_ffpredlog_ftpredlog_ttpredlog, new_p in state_ffpredlog_ftpredlog_ttpredlog_p_list:
        new_state, new_ffpredlog, new_ftpredlog, new_ttpredlog = state_ffpredlog_ftpredlog_ttpredlog
        new_fpred = kk(predlogr + new_ffpredlog)
        new_tpred = kk(predlogr + new_ftpredlog) + ((1-kk(predlogr))*kk(new_ttpredlog))
        #print('DLTIEALPHQ new_fpred={new_fpred} new_tpred={new_tpred}'.format(
        #  new_fpred=new_fpred,
        #  new_tpred=new_tpred,
        #))
        new_predlogr  = rr(ll(new_fpred/(new_fpred+new_tpred)))
        new_state_predlogr_p_list.append(
          ((new_state, new_predlogr), new_p)
        )

      state_predlogr_to_formula_data_dict[state_predlogr] = {
        'state_predlogr_p_list':new_state_predlogr_p_list,
        'else': 1
      }

      #print('RMBHBLCAAI state_predlogr={state_predlogr} new_state_predlogr_p_list={new_state_predlogr_p_list}'.format(
      #  state_predlogr=state_predlogr,
      #  new_state_predlogr_p_list=new_state_predlogr_p_list,
      #))

      for new_state_predlogr, new_p in new_state_predlogr_p_list:
        add_state_predlogr_p(new_state_predlogr, new_p*state_predlogr__p)

      # cal over
      grid_size = len(state_predlogr_to_mat_idx_dict)
      if grid_size % 10 == 0:
        a_ll_np = np.zeros((grid_size,grid_size))
        b_l_np = np.zeros((grid_size,))
        for idx0 in range(len(state_predlogr_to_mat_idx_dict)):
          a_ll_np[idx0,idx0] = 1
          state_predlogr_0 = mat_idx_to_state_predlogr_dict[idx0]
          if state_predlogr_0 not in state_predlogr_to_formula_data_dict: continue
          formula_data_0 = state_predlogr_to_formula_data_dict[state_predlogr_0]
          for state_predlogr_1, p1 in formula_data_0['state_predlogr_p_list']:
            idx1 = state_predlogr_to_mat_idx_dict[state_predlogr_1]
            a_ll_np[idx0,idx1] -= p1
          b_l_np[idx0] = formula_data_0['else']
        x_l_np = np.linalg.solve(a_ll_np, b_l_np)
        add_depth = x_l_np[0]
        if add_depth - step_limit > FLOAT_CORRECT:
          bad = True
          break

    if bad: continue
    
    grid_size = len(state_predlogr_to_mat_idx_dict)
    a_ll_np = np.zeros((grid_size,grid_size))
    b_l_np = np.zeros((grid_size,))
    for idx0 in range(len(state_predlogr_to_mat_idx_dict)):
      a_ll_np[idx0,idx0] = 1
      state_predlogr_0 = mat_idx_to_state_predlogr_dict[idx0]
      if state_predlogr_0 not in state_predlogr_to_formula_data_dict:
        a_ll_np[idx0,0] -= 1
      else:
        formula_data_0 = state_predlogr_to_formula_data_dict[state_predlogr_0]
        for state_predlogr_1, p1 in formula_data_0['state_predlogr_p_list']:
          idx1 = state_predlogr_to_mat_idx_dict[state_predlogr_1]
          #print('QKVLIRTTIW state_predlogr_1={state_predlogr_1} idx1={idx1}'.format(
          #  state_predlogr_1=state_predlogr_1,
          #  idx1=idx1,
          #))
          a_ll_np[idx0,idx1] -= p1
        b_l_np[idx0] = formula_data_0['else']

    #print(a_ll_np)
    #print(b_l_np)

    try:
      x_l_np = np.linalg.solve(a_ll_np, b_l_np)
    except:
      continue
    add_depth = x_l_np[0]
    if add_depth - step_limit > FLOAT_CORRECT:
      continue
    add_depth += 1/cross0_c_to_p_dict[c0]

    roll_data_list.append({
      'product': tg,
      'product.color': g_to_c_dict[tg],
      'parent_list': (parent_gene_0,parent_gene_1,verify_gene),
      'step_depth': 1,
      'step_count': add_depth,
      'add_depth': add_depth,
      'method': 'roll',
    })

  #print('SPEMSCKDFN roll_data_list='+str(roll_data_list))

  return roll_data_list

def cross(g0a, g0b, path_algo_cache):
  if 'cross' not in path_algo_cache:
    path_algo_cache['cross'] = {}
  ret_cache = path_algo_cache['cross']

  g00, g01 = tuple(sorted((g0a, g0b)))
  if (g00, g01) in ret_cache:
    return ret_cache[(g00, g01)]
  ret = _cross(g00, g01)
  ret_cache[(g00, g01)] = ret
  return ret

def _cross(gene0, gene1):
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


def self_cross(g, g_to_c_dict, path_algo_cache):
  if 'self_cross' not in path_algo_cache:
    path_algo_cache['self_cross'] = {}
  ret_cache = path_algo_cache['self_cross']

  if g in ret_cache:
    return ret_cache[g]
  ret = _self_cross(g, g_to_c_dict, path_algo_cache)
  ret_cache[g] = ret
  return ret

def _self_cross(g, g_to_c_dict, path_algo_cache):
  if '1' not in g:
    return g
  c = g_to_c_dict[g]
  g1_itr, _ = cross(g, g, path_algo_cache)
  g1_itr = map(lambda i:i['g'], g1_itr)
  g1_itr = filter(lambda i:i!=g, g1_itr)
  g1_itr = filter(lambda i:g_to_c_dict[i]==c, g1_itr)
  g1_itr = map(lambda i:self_cross(i, g_to_c_dict, path_algo_cache), g1_itr)
  g1_itr = list(set(g1_itr))
  if(len(g1_itr)==0): return g
  if(len(g1_itr)==1): return g1_itr[0]
  return None

def cal_merge_step(target_g, g_list, p_list, g_to_c_dict, path_algo_cache):
  target_c = g_to_c_dict[target_g]

  g0_q = collections.deque()
  for g in g_list: g0_q.append(g)
  
  g_all_set = set(g_list)
  g0_to_cs1_data_list_dict = {}
  while(len(g0_q)>0):
    g0 = g0_q.popleft()

    cs1_data_list, _ = cross(g0, g0, path_algo_cache)
    cs1_data_list = filter(lambda i:g_to_c_dict[i['g']]==target_c, cs1_data_list)
    cs1_data_list = list(cs1_data_list)
    g0_to_cs1_data_list_dict[g0] = cs1_data_list

    for cs1_data in cs1_data_list:
      g1 = cs1_data['g']
      if g1 in g_all_set: continue
      g_all_set.add(g1)
      g0_q.append(g1)

  g_all_list = list(sorted(g_all_set))
  g_to_idx_dict = { g_all_list[i]: i for i in range(len(g_all_list)) }

  p_1np = np.zeros((len(g_all_list),))
  for i in range(len(p_list)):
    g0 = g_list[i]
    p0 = p_list[i]
    g0_idx = g_to_idx_dict[g0]
    p_1np[g0_idx] = p0

  m_2np = np.zeros((len(g_all_list),len(g_all_list)))
  for g0, cs1_data_list in g0_to_cs1_data_list_dict.items():
    g0_idx = g_to_idx_dict[g0]
    for cs1_data in cs1_data_list:
      g1 = cs1_data['g']
      c1 = g_to_c_dict[g1]
      if c1 != target_c: continue
      g1_idx = g_to_idx_dict[g1]
      p = cs1_data['p']
      m_2np[g0_idx, g1_idx] = p

  target_g_idx = g_to_idx_dict[target_g]

  p_1np_sum = np.sum(p_1np)
  p_1np /= p_1np_sum
  #print('LKWFUKOTHJ p_1np={p_1np}'.format(p_1np=p_1np))

  step = 0
  while(True):
    if p_1np[target_g_idx] >= 1-ACCEPTABLE_ERR: break
    p_1np = np.matmul(p_1np, m_2np)
    p_1np_sum = np.sum(p_1np)
    p_1np /= p_1np_sum
    step += 1/p_1np_sum
    step += 1
    #print('ZJZCBKPOLK p_1np={p_1np}'.format(p_1np=p_1np))

  return step

def cross_self(g0a, g0b, g_to_c_dict, path_algo_cache):
  #print('LJFYZEYXZJ cross_self g0a={g0a} g0b={g0b}'.format(
  #  g0a=g0a,
  #  g0b=g0b
  #))

  g00, g01 = tuple(sorted((g0a, g0b)))

  cs1_data_list, _ = cross(g00, g01, path_algo_cache)

  c1_to_cs1_data_list_dict = {}
  for cs1_data in cs1_data_list:
    c1 = g_to_c_dict[cs1_data['g']]
    if c1 not in c1_to_cs1_data_list_dict:
      c1_to_cs1_data_list_dict[c1] = []
    c1_to_cs1_data_list_dict[c1].append(cs1_data)

  ret = []

  for c1, c1_cs1_data_list in c1_to_cs1_data_list_dict.items():
    if len(c1_cs1_data_list) == 1:
      cs_data = c1_cs1_data_list[0]
      ret.append({
        'product': cs_data['g'],
        'product.color': g_to_c_dict[cs_data['g']],
        'parent_list': (g00, g01),
        'add_depth': 1/cs_data['p'],
        'method': 'cross',
      })
      continue

    g1_sc_itr = c1_cs1_data_list
    g1_sc_itr = map(lambda i:i['g'],g1_sc_itr)
    g1_sc_itr = map(lambda i:self_cross(i, g_to_c_dict, path_algo_cache), g1_sc_itr)
    g1_sc_itr = set(g1_sc_itr)
    if len(g1_sc_itr) > 1: continue
    g2 = list(g1_sc_itr)[0]
    if g2 == None: continue
    if g2 == g00: continue
    if g2 == g01: continue

    target_g = g2
    g_list = list(map(lambda i:i['g'], c1_cs1_data_list))
    p_list = list(map(lambda i:i['p'], c1_cs1_data_list))

    add_step = cal_merge_step(g2, g_list, p_list, g_to_c_dict, path_algo_cache)

    ret.append({
      'product': target_g,
      'product.color': g_to_c_dict[target_g],
      'parent_list': (g00, g01),
      'add_depth': (1/sum(p_list))+add_step,
      'method': 'cross_self',
    })

  return ret

def count_one(g):
  return len(list(filter(lambda i:i=='1',g)))

def pj(j):
  return json.dumps(j, indent=2, sort_keys=True)
