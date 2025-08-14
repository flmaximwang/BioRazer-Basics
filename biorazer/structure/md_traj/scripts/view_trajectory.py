import warnings
warnings.filterwarnings('ignore')

import MDAnalysis as mda
from MDAnalysis.analysis import align
import nglview as nv
import nglview.contrib.movie as nv_movie

def view_trajectory(
        topology_path,
        trajectory_path,
        select="protein and name CA",
        ref_frame=0,
        fps=6,
        movie_output="./movie.mp4",
        movie_step=10
    ):
    '''
    Load the target MD results, align selected part to the time average structure
    Return the universe, the nglview view widget and the nglview movie object
    '''
    u = mda.Universe(topology_path, trajectory_path)
    average = align.AverageStructure(u, u, select=select, ref_frame=ref_frame).run()
    ref = average.results.universe
    aligner = align.AlignTraj(u, ref, select=select, in_memory=True).run()
    view = nv.show_mdanalysis(u)
    view._set_size("800px", "800px")
    movie = nv_movie.MovieMaker(view, output=movie_output, in_memory=True, fps=fps, step=movie_step)
    return u, view, movie