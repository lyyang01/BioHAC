
from typing import Dict, Optional, Union, List

from autogen import Agent, AssistantAgent, UserProxyAgent, MyAssistantAgent
import chainlit as cl
import os
import autogen
from autogen.coding.jupyter import LocalJupyterServer
from autogen.coding.jupyter import JupyterCodeExecutor, JupyterConnectionInfo
from autogen.agentchat import GroupChat, GroupChatManager

from get_prompt import get_prompt_plan, get_prompt_coder_hard
from chainlitagent import ChainlitAssistantAgent, ChainlitCodeAssistantAgent, ChainlitUserProxyAgent, ChainlitUserProxyAgent_new, ChainlitUserProxyAgent_skip, ChainlitAssistantAgent_ex, ChainlitUserProxyAgent_plan, ChainlitUserProxyAgent_replan
import re
import math
from typing_extensions import Annotated

#gpt config
from config import config_list



#planning
from groupchats.planning import planner, reviewer, plan_proxy, custom_speaker_selection_func_plan

#action execution
from groupchats.action_execution import coder, checker, code_proxy_hard, code_proxy_hard2, code_proxy_simple
from groupchats.action_execution import retriever_call, retriever_executor, retriever_code
from groupchats.action_execution import custom_speaker_selection_func_code

#skip
from groupchats.action_execution import skip_proxy

#decision
from groupchats.decision_maker import decisioner1, decisioner2


async def ask_helper(func, **kwargs):
    res = await func(**kwargs).send()
    while not res:
        res = await func(**kwargs).send()
    return res



import math
@cl.on_chat_start
async def on_chat_start():
    #OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')
    #try:
        #llm_config
    llm_config = {"config_list": config_list, "seed": 42, "timeout": 1000, "temperature": 0}

    ## Planner Groupchat
    #planner
    #reviewer
    #plan_proxy
    #config
    #build planner
    

    ## Skip
    #skip_proxy

    ## Refined Planner Groupchat
    #plan_refiner
    #plan_refiner_reviewer
    #user_proxy_refine

    ## Action Execution Groupchat
    #coder
    #checker
    #code_proxy_hard
    #code_proxy_hard2
    #
    #code_proxy_simple

    #========================================start chainlit server================================================================

    msg = cl.AskUserMessage(content=f"""Hello! I am a bioinformatics LLM-agent. What task would you like to get done? 
                        """, timeout=100000000, raise_on_timeout=False)
    message = await msg.send()
    
    #user input containing task and dataset
    TASK = message['output']
    
    #await cl.Message(content=f"""Starting planning on task. If you think the plan is done, please to Continue!""").send()
    
    #await cl.Message(content=f"""To plan or code?""").send()
    #await cl.make_async(decisioner1.initiate_chat)(decisioner1, message="")
    #import pdb
    #pdb.set_trace()
    #decisioner1._human_input[0] == 'plan'
    indicator_tmp = 'plan'
    #if decisioner1._human_input[0] == 'code':
    if indicator_tmp == 'code':

        PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """
        retriever_code.customized_prompt = TASK + '\n\n' + PROMPT_CODE
        code_prompt = TASK# + '\n\n' + recommend_description
        # == Begin coder ==
        indicator = 0
        MAX_ROUND = 100
        groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
        manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
        
        if len(groupchat.messages) == indicator:
            chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
            #await user_proxy.initiate_chat( manager, message=message)
        elif len(groupchat.messages) <  MAX_ROUND + indicator:
            chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
        elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
            chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
        

    else:
    
        #== Begin planner groupchat ==
        #
        #prompt = get_prompt_plan(TASK, recommend_method=recommend_description)
        prompt = get_prompt_plan(TASK, recommend_method=None) 
        #prompt = TASK + "\n\n" + "I also provide you the data details as the following.\n" +agent_data_description
        #prompt = TASK

        
        groupchat_plan = GroupChat(agents=[plan_proxy, planner, reviewer], messages=[], max_round=100, speaker_selection_method=custom_speaker_selection_func_plan)
        manager_plan = GroupChatManager(groupchat=groupchat_plan, llm_config=llm_config)

        if len(groupchat_plan.messages) == 0:
            chat_history = await cl.make_async(plan_proxy.initiate_chat)(manager_plan, message=str(prompt))
                        #await user_proxy.initiate_chat(manager, message=message)
        elif len(groupchat_plan.messages) <  100:
            chat_history = await cl.make_async(plan_proxy.send)(manager_plan, message=TASK,)
        elif len(groupchat_plan.messages) ==  100:
            chat_history = await cl.make_async(plan_proxy.send)(manager_plan, message="exit",)
        
        #== End planner groupchat ==

        final_plan = groupchat_plan.messages[-1]['content']
        plan_str = final_plan

        await cl.Message(content=f"""To sub-steps or code?""").send()
        await cl.make_async(decisioner2.initiate_chat)(decisioner2, message="")
        #import pdb
        #pdb.set_trace()

        if decisioner2._human_input[0] == 'code':
            PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """
            retriever_code.customized_prompt = TASK + '\n\n' + PROMPT_CODE
            #goal_tmp = "write codes for the whole plan, not just one step. If necessary, please use GPU settings in your code to accelerate the training speed."
            code_prompt = plan_str#get_prompt_coder_hard(plan_str, goal_tmp, None)
            # == Begin coder ==
            indicator = 0
            MAX_ROUND = 100
            groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
            manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
            
            if len(groupchat.messages) == indicator:
                chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
                #await user_proxy.initiate_chat(manager, message=message)
            elif len(groupchat.messages) <  MAX_ROUND + indicator:
                chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
            elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
                chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
        
        else:
            #import pdb
            #pdb.set_trace()
            sub_steps = plan_str.split('\n\n')#[1:-1]
            sub_steps_temp = []
            for i in sub_steps:
                pattern = r'Step \d'
                #if 'Step' in i:
                if re.search(pattern, i):
                    sub_steps_temp.append(i)
            sub_steps = sub_steps_temp

            
            indicator = 0
            previous_step = ""
            coder_messages = []
            tmp_coder_messages = ""
            for idx, sub_step in enumerate(sub_steps):
                ori_sub_step = sub_step

                #=====start skip===
                message2 = "Step {} is:\n".format(str(idx+1)) + sub_step + "\n\nYour goal: determine whether skip the current step or refined planning."
                await cl.Message(content=f"""Starting Skipping. Please decide whether to skip!""").send()
                await cl.make_async(skip_proxy.initiate_chat)(skip_proxy, message=message2)
                # == End skip ==
                
                indicator_ = skip_proxy._human_input

                if indicator_[0] == 'skip-step':
                    continue

                goals = "write the code for Step {}.".format(str(idx+1))
                #message2 ="the user task need: " + TASK + '.\n\n' + "The total plan:\n" + plan_str + "\n\n" + "The step {} is:\n".format(str(idx+1)) + sub_step# + "\n\n" + "Your goal: write the code for the step {}.'".format(str(idx+1))
                if previous_step == "":
                    message2 ="Step {} is:\n".format(str(idx+1)) + sub_step# + "\n\n" + "Your goal: write the code for the step {}.'".format(str(idx+1))
                else:
                    #message2 ="Previous steps are:\n"+ previous_step + "\n" + "Codes of previous steps are:\n"  + tmp  + "\n\n" +"Step {} is:\n".format(str(idx+1)) + sub_step
                    message2 ="Codes of previous steps are:\n"  + tmp  + "\n\n" +"Step {} is:\n".format(str(idx+1)) + sub_step

                code_prompt = get_prompt_coder_hard(message2, goals, None)
                
                
                PROMPT_CODE = """Now, you're a retrieve augmented Coder. You goal is to write the code for the current step task based on your own knowledge and the context. Context generally contains how specific packages are used. Thus, you need to extract useful information from the context.\nNote that you do not need to follow the code in Context strictly. In other words, You can modify or discard the lines of code in the context that are not necessary for the task, or you can write new code yourself.\nIf Context does not contain the knowledge you need to write the code, then ignore Context and write the code using only your own knowledge.\n\nContext is: {input_context}
                """

                retriever_code.customized_prompt = "the task of step {}:\n".format(str(idx+1)) + ori_sub_step + '\n\n' + PROMPT_CODE

                # == Begin coder ==
                MAX_ROUND = 100
                groupchat = GroupChat(agents=[code_proxy_hard, code_proxy_hard2, coder, checker, retriever_executor], messages=[], max_round=MAX_ROUND, speaker_selection_method=custom_speaker_selection_func_code)
                manager = GroupChatManager(groupchat=groupchat, llm_config=llm_config)
                
                if len(groupchat.messages) == indicator:
                    chat_history = await cl.make_async(code_proxy_hard2.initiate_chat)(manager, message=str(code_prompt)+ '\n\nretrieve content')
                    #await user_proxy.initiate_chat( manager, message=message)
                elif len(groupchat.messages) <  MAX_ROUND + indicator:
                    chat_history = await cl.make_async(code_proxy_hard.send)(manager, message=TASK,)
                elif len(groupchat.messages) ==  MAX_ROUND + indicator:  
                    chat_history = await cl.make_async(code_proxy_hard.send)(manager, message="exit",)
                # == End code ==
                
                #indicator = len(groupchat.messages)
                
                previous_step = previous_step + sub_step + "\n"

                
                
                for i in groupchat.messages:
                    if i['name'] == "Coder":
                        if "```python" in i['content']:
                            pattern = r"```python(.*?)```"
                            matches = re.findall(pattern, i['content'], re.DOTALL)
                            #import pdb
                            #pdb.set_trace()
                            #code_tmp = i["content"].split()
                            
                            #tmp_coder_messages = i["content"]
                            tmp_coder_messages = '```python' + '\n'.join(matches) + '```' + '\n' 

                            
                coder_messages.append(tmp_coder_messages)
                tmp = '\n'.join(coder_messages)
        
        
    f = open('./chat_result/multi-level-clustering.txt', 'w')
    f.write(str(groupchat.messages))
    f.close()
