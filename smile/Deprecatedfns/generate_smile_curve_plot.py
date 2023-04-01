def generate_smile_curve_plot(sa_deg):
    plt.figure(figsize=(50, 20))
    for i in range(np.shape(sa_deg)[0]):
      plt.plot(range(10), sa_deg[i], marker='+',linestyle='-');

    plt.show()