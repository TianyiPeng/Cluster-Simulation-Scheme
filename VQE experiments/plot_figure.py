
# We import plotting tools 
import matplotlib.pyplot as plt 
from   matplotlib import cm
from   matplotlib.ticker import LinearLocator, FormatStrFormatter
import pickle

fw = open('SPSA_result_prepared_plot_200.txt', 'rb')
[e_minus, e_plus, e_ideal, e_nature, e_exact] = pickle.load(fw)
print(len(e_minus))

T = range(len(e_minus))

plt.plot(T, e_ideal, 'g-', T, e_nature, 'm--', T, e_exact, 'k--')
leg = plt.legend([ 
            r'performance of $\theta(k)$', 'Kandala et al. [44]', 'exact energy'])
leg.get_frame().set_linewidth(1.5)
leg.get_frame().set_edgecolor('k')
ax = plt.gca()
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)

plt.xlabel('Iteration, $k$')
plt.ylabel('Energy (hartree)')
plt.savefig('VQE_no_SPSA.eps', format='eps')
plt.show()
