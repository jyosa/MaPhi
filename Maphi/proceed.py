import os

'''
    check if exist errors in molecules
'''    


def ask_user():
    check = str(input("There are some molecules with problems, are you going to fix it? or can we cotinue with the calculation?\n\
type  (Y/N): ")).lower().strip()
    try:
        if check[0] == 'y':
            return True
        elif check[0] == 'n':
            return False
        else:
            print('Invalid Input')
            return ask_user()
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return ask_user()


def check_err(molepath):
    dir_n = str(molepath + "/" + "molecules_error")
    if os.path.exists(dir_n) and os.path.isdir(dir_n):
        if not os.listdir(dir_n):
            print("\n\n checking for errors...: All molecules are OK!")
        else:    
            ans = ask_user()
            if ans == True:
                print("\n\n Cool! We are going to proceed.\n")
            else:
                print("\nFix molecules and run again Maphi!!")
                exit(1)            
 
 