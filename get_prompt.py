def get_prompt_plan(user_input, agent_data_description=None, recommend_method=None):
        prompt = {
                    "role": "Act as a bioinformatician, following these rules strictly!",
                    "rules": {
                        "Assume all necessary packages are installed; do not include installation steps.",
                        "If new data is generated, save it as a new file.",
                        "When modifying the plan, keep previous details such as data paths and variable names.",
                        "Avoid including specific code blocks in your plan.",
                        "After receiving suggestions from 'Reviewer', make changes directly to the previous plan.",
                        "Start the planning structure with 'Plan for ...', presenting steps as 'Step 1... Step 2...'.",
                        "Provide specific file paths in each plan step."
                    },
                    "input": [
                        "I provide all user input, which may include task and data information.",
                        user_input
                    ],
                    "recommend_method": [
                        "I have provided recommended methods. Choose one for the plan; if none are available, use your own knowledge."
                    ]
                }

        return prompt

def get_prompt_coder_hard(msg, goal, data_list):
        prompt = {
                    "role": "Act as a bioinformatician, the rules must be strictly followed!",
                    "rules": [
                        #"All rules must be followed strictly!",
                        #"You should write code for all steps and do not omit the code like # Step 2: xxx...!",
                        #"All code should be put in one code block.",
                        #"You should put code in different code blocks basd on the steps in the plan.",
                        "Assume that the environment configuration is all ready",
                        "You must give answer with code blocks. Specifically, your answer should include '```python' and '```'",
                        "Assuming you have a jupyter kernel, so you can continuously execute code cell by cell.",
                        "You should write code for the current one step, assuming that the code for the previous steps has already been written.",
                        #"When you need to write a new file, you can first determine whether the file exists on the front end of the code. If it exists, the code execution of the current step can be skipped. You can use the if else statement to complete this part.",                       
                        #"Include necessary visualizations in the code!",
                        "You need to use Python to write your code, and for some terminal commands, please use 'subprocess' to implement them in Python",
                        "For functions used in your code, use your own knowledge to add essential parameters for the function. You can also use your own knowledge to fill in missing code steps or ignore unnecessary code steps",
                        #"In addition to the information you get from the 'plan_refiner' or 'Planner', you can also use your own knowledge to fill in missing code steps or ignore unnecessary code steps.",
                        #"You need to wrap the code block in ```python and ```",
                    ],
                    "input": ["I provide all information to you including codes of previous steps and the current step description", msg],
                    #"data_information": [
                    #        "I also provide the data path and data description to you.",
                    #        data_list
                    #    ],
                    "goal": ["I provide your goal", goal],
            }
        return prompt
