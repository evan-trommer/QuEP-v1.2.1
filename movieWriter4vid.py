import os
import pathlib
import cv2
import numpy as np
from PIL import Image
import glob
import time
 
def generatemovie4vid(fps,image_folder_1,image_folder_2,image_folder_3,image_folder_4):

    img_array = []
    #path = os.path.join(r'/gpfs/home/nmanzella/PWA/QuEP')

    
    dim_scale = 0.5
    

    #fps = 8
    print(f"Running at {fps} frames per second")

    
    dim_scale_ = "{:03.2f}".format(dim_scale).replace(".","-")
    video_name = f'00-progression-movie-4vid__{dim_scale_}res__{fps}fps.mp4'
    
    # GET FIRST IMAGES
    os.chdir(image_folder_1)

    images1 = [img1 for img1 in os.listdir(image_folder_1)
            if img1.endswith(".jpg") or
                img1.endswith(".jpeg") or
                img1.endswith("png")]

    images1.sort()# = sorted(images1)
    
    # GET SECOND IMAGES
    os.chdir(image_folder_2)

    images2 = [img2 for img2 in os.listdir(image_folder_2)
            if img2.endswith(".jpg") or
                img2.endswith(".jpeg") or
                img2.endswith("png")]

    images2.sort()# = sorted(images2)

    # GET THIRD IMAGES
    os.chdir(image_folder_3)

    images3 = [img3 for img3 in os.listdir(image_folder_3)
            if img3.endswith(".jpg") or
                img3.endswith(".jpeg") or
                img3.endswith("png")]

    images3.sort()# = sorted(images3)

    # GET FOURTH IMAGES
    os.chdir(image_folder_4)

    images4 = [img4 for img4 in os.listdir(image_folder_4)
            if img4.endswith(".jpg") or
                img4.endswith(".jpeg") or
                img4.endswith("png")]

    images4.sort()# = sorted(images4)

    print(f"\nUsing {len(images1)} images from: {image_folder_1}")
    print(f"\nUsing {len(images2)} images from: {image_folder_2}")
    print(f"\nUsing {len(images3)} images from: {image_folder_3}")
    print(f"\nUsing {len(images4)} images from: {image_folder_4}")

    # Array images should only consider
    # the image files ignoring others if any
    filenumber = "{:05.1f}".format(0).replace(".","-") 
    img1 = cv2.imread(os.path.join(image_folder_1, f'progression-x-{filenumber}mm.png'))
    img2 = cv2.imread(os.path.join(image_folder_2, f'progression-x-{filenumber}mm.png'))
    cropped_img1 = img1#[:2415,:]
    cropped_img2 = img2#[250:,:]

    img3 = cv2.imread(os.path.join(image_folder_3, f'progression-x-{filenumber}mm.png'))
    img4 = cv2.imread(os.path.join(image_folder_4, f'progression-x-{filenumber}mm.png'))
    cropped_img3 = img3#[:2415,:]
    cropped_img4 = img4#[250:,:]

    vis_left = np.concatenate((cropped_img1, cropped_img2), axis=0)
    vis_right = np.concatenate((cropped_img3, cropped_img4), axis=0)
    vis = np.concatenate((vis_left, vis_right), axis=1)

    old_height, old_width, layers = vis.shape
    new_width = int(old_width * dim_scale)
    new_height = int(old_height * dim_scale)
    print(f"New dimensions for video: {new_width}x{new_height}px")
    frame = cv2.resize(vis, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)

    

    fourcc = cv2.VideoWriter_fourcc(*"h264")
    #0x00000021 #Seawulf: 0x31637661
    video = cv2.VideoWriter(video_name, 0x00000021, fps, (new_width, new_height)) 

    start_time_video = time.time()
    t = time.localtime()
    curr_time = time.strftime("%H:%M:%S", t)
    print("\nMovie generation - START TIME: ", curr_time)

    # Appending the images to the video one by one
    for i in range(0,len(images1)):
        filenumber = "{:05.1f}".format(i).replace(".","-") 
        img1 = cv2.imread(os.path.join(image_folder_1, f'progression-x-{filenumber}mm.png'))
        img2 = cv2.imread(os.path.join(image_folder_2, f'progression-x-{filenumber}mm.png'))
        cropped_img1[i] = img1#[:2415,:]
        cropped_img2[i] = img2#[250:,:]

        img3 = cv2.imread(os.path.join(image_folder_3, f'progression-x-{filenumber}mm.png'))
        img4 = cv2.imread(os.path.join(image_folder_4, f'progression-x-{filenumber}mm.png'))
        cropped_img3[i] = img3#[:2415,:]
        cropped_img4[i] = img4#[250:,:]

        vis_left = np.concatenate((cropped_img1, cropped_img2), axis=0)
        vis_right = np.concatenate((cropped_img3, cropped_img4), axis=0)
        vis = np.concatenate((vis_left, vis_right), axis=1)
        frame = cv2.resize(vis, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)

        video.write(frame) 

    # Deallocating memories taken for window creation
    cv2.destroyAllWindows() 
    video.release()  # releasing the video generated


    print(f"Movie generated! File stored at {os.path.join(image_folder_2,video_name)}")
    t_video_end = time.localtime()
    curr_time_video_end = time.strftime("%H:%M:%S", t_video_end)
    print("Movie generation - END TIME: ", curr_time_video_end)
    print("Movie generation - END TIME: ", (time.time() - start_time_video)/60, " min\n")

