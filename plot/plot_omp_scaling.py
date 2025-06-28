import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('omp_scaling_results/omp_thread_times.csv', header=None)
df.columns = ['Threads', 'Min', 'Max', 'Avg']

# Plot Min, Max, Avg Thread Computation Time
plt.figure(figsize=(10, 6))
plt.plot(df['Threads'], df['Min'], marker='o', label='Min Thread Time')
plt.plot(df['Threads'], df['Avg'], marker='o', label='Avg Thread Time')
plt.plot(df['Threads'], df['Max'], marker='o', label='Max Thread Time')
plt.xlabel('Number of Threads')
plt.ylabel('Computation Time (s)')
plt.title('Thread Computation Time vs Number of Threads')
plt.xscale('log', base=2)
plt.xticks(df['Threads'])
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.savefig('thread_times.png', dpi=300)
plt.show()