# permet de définir les couleurs d'affichage
set palette defined ( 0 '#F7FBFF',\
                      1 '#DEEBF7',\
                      2 '#C6DBEF',\
                      3 '#9ECAE1',\
                      4 '#6BAED6',\
                      5 '#4292C6',\
                      6 '#2171B5',\
                      7 '#084594' )

# pour faire des png. Commenter pour un affichage classique.
set terminal png

# forcer l'étendue de la colorbar. A changer selon les cas.
set cbrange [-0.2:1.2]

# force la mise à l'échelle des axes.
set size ratio -1

# si affichage 3D, force l'étendue de l'axe des z. A changer selon les cas.
set zrange [-0.2:1.2]
# si affichage 3D, fixe le point de vue. A changer selon les cas.
set view 48,132

do for [i = 0:250] {
    t=i*0.01
    # Set output file name
    set output "sol.".sprintf("%05d", i).".png"
    # Set title
    set title "t = ".sprintf("%f", t)." s"." (i = ".sprintf("%d", i).")"
    # Display title
    show title

    # Configuration for a 3D surface plot
    set pm3d interpolate 2,2
    set dgrid3d 50,50,2
    set view 60, 30, 1, 1

    # Plotting a 3D surface
    splot "./sol.".i."_.dat" using 1:2:3 with pm3d
    

    # Uncomment the following line to pause between frames (useful for animation viewing)
}
