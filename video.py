import cv2
import os

# Chemin vers le dossier contenant les images
folder_path = '/home/yassir/Desktop/results/chp_projet'
# Nom du fichier vidéo de sortie
video_name = 'output_video.mp4'

images = [img for img in os.listdir(folder_path) if img.endswith(".png")]
images.sort()  # Assurez-vous que les images sont triées, si nécessaire

# Capturez les dimensions d'une image
frame = cv2.imread(os.path.join(folder_path, images[0]))
height, width, layers = frame.shape

# Création de l'objet vidéo
video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'mp4v'), 10, (width, height))

for image in images:
    video.write(cv2.imread(os.path.join(folder_path, image)))

cv2.destroyAllWindows()
video.release()
