# index-mp.py can be used to propagate objects after simulation through fields in main.py
# Available features include weighting of particles and masking of regions

# Include file imports
from logging import raiseExceptions
import sys
import time
import importlib
import numpy as np
import include.showQuickEvolution as showEvol_Q
import include.showFullEvolution as showEvol_F
import include.makeFullAnimation as makeFullAni
import include.viewProbe as viewProbe
import include.writeFullEvolData as writeHist
import include.weighting_masks_function_rprism as weightmaskFunc
import include.plotWeights as plotWeights
import include.findFocalY as findFocalY
import include.plot2DTracks as plot2D
import include.plot3DTracks as plot3D
import include.findWaist as findWaist
import multiprocessing as mp
import include.movieWriter as movieWriter
#import include.movieWriter4vid as movieWriter4v
import tqdm
import pickle
from DebugObjectModule import DebugObject
from random import randint



# Be sure to change .npz file name location from main.py output!
# Put .npz file in /data directory

# Weighting Options (Only applicable for showFullEvolution and makeFullAnimation plot):
useWeights_x = False                 # NOT CURRENTLY IN USE - LEAVE FALSE - Use weights in x-direction
useWeights_y = False                  # Use gaussian weights in y-direction
useWeights_xi = False                 # Use gaussian weights in xi-direction

skipWeightingCalc = False            # Skip weighting calculation and use imported pre-calculated weights
saveWeights = False                 # Save weights to .npz file (Remember to move to ./data directory!)

# Masking Options:
useMasks_xi = False                 # Use masks in xi-direction (Vertical; done during weighting)
useMasks_y = False                  # Use masks in y-direction (Horizontal; done during weighting)
useMasks_x = False                  # NOT CURRENTLY IN USE - LEAVE FALSE - Use masks in x-direction (transverse; done during weighting)

# Plotting Scripts
showQuickEvolution = False           # View evolution of probe after leaving plasma at inputted x_s in scatter plots # Use for low density probes
showFullEvolution = False             # View full evolution of probe at hardcoded locations in colored histograms # Use for high density probes
makeFullAnimation = False#True
writeHistData = False

# Gaussian Weighting Testing
plotWeightsx = False                  # Plot w vs xi (ONLY for single line of particles in x-dir)
plotWeightsy = False                  # Plot w vs y (ONLY for single line of particles in y-dir)
plotWeightsxi = False                  # Plot w vs y (ONLY for single line of particles in xi-dir)
plotWeightsxiy = False                 # Plots initial particle density map
plotWeights3D = False                  # Plots y, xi, and 2D cross section

# DEBUG PLOTTING
plot2DTracks = True#False                 # View 2D projections of trajectories (SET ALL OTHERS TO FALSE & ONLY USE FOR SINGLE PARTICLE)
findFocal = False
plot3DTracks = False
findW = False

#moreMovies = True

if __name__ == '__main__':
    # Start of main()
    # Initialize multiprocessing.Pool()
    numberOfCores = 6 #mp.cpu_count() 
    print(f"Number of cores used for multiprocessing: {numberOfCores}")
    pool = mp.get_context('spawn').Pool(numberOfCores)
    if (len(sys.argv) >= 2):
        
        # Begin timing index file runtime
        start_time = time.time()
        t = time.localtime()
        curr_time = time.strftime("%H:%M:%S", t)
        print("index.py - START TIME: ", curr_time)
        
        # Get inital conditions of probe again
        input_fname_1 = str(sys.argv[1])
        print("Using initial conditions from ", input_fname_1)
        init = importlib.import_module(input_fname_1)
        sim_name = init.simulation_name
        shape_name = init.shape
        xden = init.xdensity
        yden = init.ydensity
        xiden = init.xidensity
        res = init.resolution
        iter = init.iterations
        mode = init.mode
        fname = init.fname
        debugmode = init.debugmode
        x_c = init.x_c
        y_c = init.y_c
        xi_c = init.xi_c
        px_0 = init.px_0
        py_0 = init.py_0
        pz_0 = init.pz_0
        x_s = init.x_s
        s1 = init.s1
        s2 = init.s2

        if len(sys.argv) == 3:
            # Get initial conditions of beam
            input_fname_2 = str(sys.argv[2])
            print("Using beam conditions from ", input_fname_2)
            beaminit = importlib.import_module(input_fname_2)
            beamx_c = beaminit.beamx_c
            beamy_c = beaminit.beamy_c
            beamxi_c = beaminit.beamxi_c
            sigma_x=beaminit.sigma_x
            sigma_y=beaminit.sigma_y
            sigma_xi=beaminit.sigma_xi
        else:
            print("WARNING: No gaussian weights inputted. Make sure not using weighting!")

        # Load data from npz file export from main.py
        data = np.load('./data/' + fname) # Change this line as needed
        x_0 = data['x_init']
        y_0 = data['y_init']
        xi_0 = data['xi_init']
        z_0 = data['z_init']
        x_f = data['x_dat']
        y_f = data['y_dat']
        xi_f = data['xi_dat']
        z_f = data['z_dat']
        px_f = data['px_dat']
        py_f = data['py_dat']
        pz_f = data['pz_dat']
        t0 = data['t_dat']
        
        if debugmode == True:
            # Load debug data from .obj file export from main.py
            file = open("./data/"+fname[:-4]+"-DEBUG.obj", 'rb') 
            debug = pickle.load(file)[0]
            file.close
            print(debug)
            print(type(debug))
            #print(debug.x_dat)
            x_dat = debug.x_dat
            y_dat = debug.y_dat
            z_dat = debug.z_dat
            xi_dat = debug.xi_dat
            Fx_dat = debug.Fx_dat
            Fy_dat = debug.Fy_dat
            Fz_dat = debug.Fz_dat
            px_dat = debug.px_dat
            py_dat = debug.py_dat
            pz_dat = debug.pz_dat
            
        noObj = len(x_0) # Number of particles in the simulation (2D Projection)

        # WEIGHTING IMPORTS/SAVING
        rand = "{:02d}".format(randint(0,99))
        weights_fname = fname[:-4] + "-weights-" + rand
        if (skipWeightingCalc):
            print("Not using weighting!")
            #Do not use these when simulating probe with no weights - i.e. for plot2Dtracks
            data = np.load('./data/' + weights_fname + '.npz') # Change this line as needed
            w = data['w']
            print(f"\nUsing weights from {'./data/' + weights_fname + '.npz'}...\n")
        else:
            # Create weighting array with appropriate masks
            w = []
            w = [1 for k in range(0,noObj)] #Creates default array of weights with length noObj, each with value 1
                  
            start_time_w = time.time()
            t_w = time.localtime()
            curr_time_w = time.strftime("%H:%M:%S", t_w)
            print("\nWeighting calculations - START TIME: ", curr_time_w)

            # Call weighting function getWeights 
            # Note: w_virt, xv, yv, xiv, only used for debugging purposes
            #w, w_export1, w_y, w_xi = weightmaskFunc.getWeights(beamx_c,beamy_c,beamxi_c,x_c,y_c,xi_c,s1,s2,xden,yden,xiden,res,sigma_x,sigma_y,sigma_xi,noObj,t0,useWeights_x,useWeights_y,useWeights_xi,useMasks_x,useMasks_xi,useMasks_y)   
            
            
            t_w_end = time.localtime()
            curr_time_w_end = time.strftime("%H:%M:%S", t_w_end)
            print("Weighting calculations - END TIME: ", curr_time_w_end)
            print("Weighting calculations - DURATION: ", (time.time() - start_time_w)/60, " min\n")
            
            if (saveWeights):
                np.savez(weights_fname, w=w)
                print(f"\nWeights saved to {weights_fname + '.npz'}\n") #Saves weights for reuse

        # Plot data points
        print("Plotting...")
        if (showQuickEvolution):
            showEvol_Q.plot(x_f, y_f, xi_f, z_f, px_f, py_f, pz_f, sim_name, shape_name, x_s, noObj, iter) # Note: does not use weights
        
        if (showFullEvolution):
            showEvol_F.plot(x_f, y_f, xi_f, z_f, px_f, py_f, pz_f, w, sim_name, shape_name, noObj, iter)
        
        if (makeFullAnimation):
            #Prepare plotting variables
            plasma_bnds, slices, xs_norm, yslice, zslice, bin_edges_z, bin_edges_y, cmap, cmin, vmin_, vmax_, zmin, zmax, ymin, ymax, fps, new_path, screen_dists = makeFullAni.prepare(sim_name, shape_name, noObj, rand)
            
            # Multiprocessing: propagate to each screen and create frame
            start_time_pfc = time.time()
            t_pfc = time.localtime()
            curr_time_pfc = time.strftime("%H:%M:%S", t_pfc)
            print("Multiprocessing propagation and frame creation - START TIME: ", curr_time_pfc)
            
            pool.starmap(makeFullAni.plotmp,[(i,x_f,y_f,z_f,px_f,py_f,pz_f, w, xden, plasma_bnds, xs_norm, yslice, zslice, bin_edges_z, bin_edges_y, cmap, cmin, vmin_, vmax_, zmin, zmax, ymin, ymax, new_path, screen_dists) for i in range(0,slices)])

            pool.close()

            pool.join()    
                      
            t_pfc_end = time.localtime()
            curr_time_pfc_end = time.strftime("%H:%M:%S", t_pfc_end)
            print("MP PFC - END TIME: ", curr_time_pfc_end)
            print("MP PFC - DURATION: ", (time.time() - start_time_pfc)/60, " min\n")
            
            # Stitch frames into movie
            #movieWriter.generatemovie(fps,new_path)
            
        if (writeHistData):
            writeHist.plot(x_f, y_f, xi_f, z_f, px_f, py_f, pz_f, sim_name, shape_name, noObj, iter)
        
        if (plotWeightsy):
            plotWeights.ploty(w_y, x_0, y_0, xi_0, z_0, s1, s2, yden, xiden, beamy_c, sigma_y)
        
        if (plotWeightsx):
            #plotWeights.plotx(w, x_0, y_0, xi_0, z_0, s1, s2, beamx_c,beamy_c,beamxi_c,sigma_x,sigma_y,sigma_xi)
            raise NotImplementedError("This functionality is not currently implemented")
            # This function could be easily added when needed

        if (plotWeightsxi):
            plotWeights.plotxi(w_xi, x_0, y_0, xi_0, z_0, s1, s2, yden, xiden, beamxi_c, sigma_xi)

        if (plotWeightsxiy):
            # Plots initial map of particle density
            plotWeights.plotweightsxiy(y_0,xi_0, w, rand)
                
        if (plotWeights3D):
            plotWeights.plotcross(w_export1, x_0, y_0, xi_0, z_0, s1, s2, yden, xiden, beamxi_c, sigma_xi)
            #plotWeights.ploty(w_y, x_0, y_0, xi_0, z_0, s1, s2, yden, xiden)
            #plotWeights.plotxi(w_xi, x_0, y_0, xi_0, z_0, s1, s2, yden, xiden)

        if (findFocal):
            findFocalY.calculate(x_0, y_0, xi_0, z_0, x_dat, y_dat, z_dat, xi_dat, px_f, py_f, pz_f, sim_name, shape_name, x_s, s1, s2)
        
        if (plot2DTracks):
            print("Plotting 2D Tracks...")          
            plot2D.plot(x_dat, y_dat, z_dat, xi_dat, Fx_dat, Fy_dat, Fz_dat, px_dat, py_dat,pz_dat, sim_name, shape_name, s1, s2, noObj, fname)           
                      
        if (plot3DTracks):
            plot3D.plot(x_dat,y_dat,z_dat,xi_dat,sim_name,shape_name,s1,s2,noObj)
        
        if (findW):
            findWaist.calculate(x_0,y_0,xi_0,z_0,x_dat,y_dat,z_dat,xi_dat,px_f,py_f,pz_f,sim_name,shape_name,x_s,s1,s2)
     
                      
        # End timing index file runtime
        tf = time.localtime()
        curr_time_f = time.strftime("%H:%M:%S", tf)
        print("index.py - END TIME: ", curr_time_f)
        print("index.py - DURATION: ", (time.time() - start_time)/60, " min\n")
    else:
        print("Improper number of arguments. Expected 'python3 index.py <fname> *<fname>'")
        exit()                 
    
    pool.close()

    pool.join()