import matplotlib.pyplot as plt

def SMM_bar_graph(trs,ref,abs):
    data = {'Transmission':trs, 'Reflection':ref, 'Absorption':abs}
    smm_results = list(data.keys())
    values = list(data.values())
    plt.bar(smm_results,values,width=0.4)
    plt.ylabel('% Power')
    plt.ylim(0,100)