
def selfies_split(selfies):
    return selfies.replace(']', '] ').split()

def find_int(token):
    for character in token:
                if character.isdigit():
                    return int(character)
                else:
                    pass

def get_numeric_value(token):
    num_dict={
        "[C]" : 0,
        "[Ring1]" : 1,
        "[Ring2]" : 2,
        "[Branch1]" : 3,
        "[=Branch1]" : 4,
        "[#Branch1]" : 5,
        "[Branch2]" : 6,
        "[=Branch2]" : 7,
        "[#Branch2]" : 8,
        "[O]" : 9,
        "[N]" : 10,
        "[=N]" : 11,
        "[=C]" : 12,
        "[#C]" : 13,
        "[S]" : 14,
        "[P]" : 15
    }
    value = num_dict[token]
    return value

def calculate_len_connect(i,nb_num_token,selfies_token):
    num=0
    skip_idx=[]
    for j in range(1,nb_num_token+1):
        skip_idx.append(i+j)
        num+=get_numeric_value(selfies_token[i+j]) * (16**(nb_num_token-j))
    num+=1
    return num,skip_idx

def token_classify (selfies_token):
    
    ### initialization of lists
    atoms_idx_list=[]
    branch_idx_list=[]
    ring_idx_list=[]
    skip_list=[]

    for i,token in enumerate(selfies_token):
        if i in skip_list:
            continue

        if token.find("Branch") !=-1:
            nb_num_token=find_int(token)
            num,skip_idx=calculate_len_connect(i,nb_num_token,selfies_token)
            skip_list.extend(skip_idx)
            branch_idx_list.append([i,nb_num_token,num])

        elif token.find("Ring")!=-1:
            nb_num_token= find_int(token)
            num,skip_idx=calculate_len_connect(i,nb_num_token,selfies_token)
            skip_list.extend(skip_idx)
            ring_idx_list.append([i,nb_num_token,num])

        else:
          atoms_idx_list.append(i)  
    
    return atoms_idx_list,branch_idx_list,ring_idx_list

def find_start_ring(start_token_idx,ring_length,atoms_idx_list):
    start_position=-1
    if start_token_idx> atoms_idx_list[-1]:
        start_index=len(atoms_idx_list)-1-ring_length-1
        start_position=atoms_idx_list[start_index]
        return start_position
    else:
        for relative_index,atoms_index in enumerate(atoms_idx_list):
            if atoms_index< start_token_idx:
                pass
            else:
                start_index=relative_index-ring_length-1
                start_position=atoms_idx_list[start_index]
                break
        return start_position

def find_ring_coord(ring_idx_list,atoms_idx_list):
    pre_ring_block=[]
    for ring_group in ring_idx_list:
        start_token_idx=ring_group[0]
        nb_num_token=ring_group[1]
        ring_length=ring_group[2]
        start_position=find_start_ring(start_token_idx,ring_length,atoms_idx_list)
        end_position=start_token_idx +nb_num_token
        pre_ring_block.append([start_position,end_position])
    return pre_ring_block 

def find_branch_coord(branch_idx_list):

    pre_branch_block=[]
    for branch_group in branch_idx_list:
        
        start_token_idx=branch_group[0]
        nb_num_token=branch_group[1]
        branch_length=branch_group[2]
        
        end_position=start_token_idx+nb_num_token+ branch_length
        start_position=start_token_idx-1

        if pre_branch_block:
            if pre_branch_block[-1][1] == start_position:
                start_position=pre_branch_block[-1][0]
                pre_branch_block.pop()
            else:
                pass
        else:
            pass
        pre_branch_block.append([start_position,end_position])
    
    return pre_branch_block  

class Solution: 
   def solve(self, intervals):            
       intervals.sort() 
       ans = [] 
       for s, e in intervals: 
         if ans and s <= ans[-1][1]: 
            ans[-1][1] = max(ans[-1][1], e) 
         else: 
            ans.append([s, e]) 
       return ans 

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

def get_index_main_chain(len_selfies,pre_block_list):
  idx_all_selfies=list(range(len_selfies))
  all_blocks_idx=[]
  for block in pre_block_list:
    temp=list(range(block[0],block[1]+1) ) 
    all_blocks_idx.extend(temp)
  # print(all_blocks_idx)
  all_blocks_idx=set(all_blocks_idx)
  main_chain_idx=[i for i in idx_all_selfies if i not in all_blocks_idx]
  # print(main_chain_idx)
  if  len(main_chain_idx)==0 :
    return None
  elif len(main_chain_idx)==1:
    main_chain_block=[[main_chain_idx[0],main_chain_idx[0]]]
    return main_chain_block

  else:
    temp=list(group(main_chain_idx))
    main_chain_block= [list(tuple_) for tuple_ in temp]
    # print(main_chain_block)
    return main_chain_block

def find_block(selfies_token,atoms_idx_list,branch_idx_list,ring_idx_list):
  pre_branch_block=find_branch_coord(branch_idx_list)
  # print(f'branch block: {pre_branch_block} ')
  pre_ring_block=find_ring_coord(ring_idx_list,atoms_idx_list)
  # print(f'ring block: {pre_ring_block} ')
  pre_block= pre_branch_block + pre_ring_block
  ob= Solution()
  pre_block_list=ob.solve(pre_block)
  # print(f'all blocks: {pre_block_list} ')
  main_blocks=get_index_main_chain(len(selfies_token),pre_block_list)
  if main_blocks :
    # print(main_blocks)
    block_list=pre_block_list +  main_blocks
    orderded_block_list=ob.solve(block_list)
    selfies_frag_list=["".join(selfies_token[block[0]:block[1]+1]) for block in orderded_block_list]
    return selfies_frag_list
  elif len(pre_block)>1:
    orderded_block_list=pre_block_list
    selfies_frag_list=["".join(selfies_token[block[0]:block[1]+1]) for block in orderded_block_list]
    return selfies_frag_list    
  else:
    return None

def process_selfies(selfies):
    tokens=selfies_split(selfies)
    # print(tokens)
    atoms_idx_list,branch_idx_list,ring_idx_list=token_classify(selfies_token=tokens)
    # print(f'branch: {branch_idx_list}')
    # print(f'ring: {ring_idx_list}')
    selfies_frag=find_block(tokens,atoms_idx_list,branch_idx_list,ring_idx_list)
    return selfies_frag
