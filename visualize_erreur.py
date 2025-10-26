import matplotlib.pyplot as plt
import numpy as np

# Vos données
x = np.array([0.2, 0.10, 0.05])
y = np.array([0.13, 0.078, 0.069])

# Calcul du logarithme de x et y
log_x = np.log(x)
log_y = np.log(y)

# Création du graphique
plt.figure(figsize=(8, 6))
plt.plot(log_x, log_y, marker='o')  # Tracer les points avec des marqueurs
plt.xlabel('log(x)')
plt.ylabel('log(y)')
plt.title('Log-log plot of y versus x')
plt.grid(True)  # Ajouter une grille pour faciliter la lecture
plt.show()
