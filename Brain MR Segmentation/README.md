# FSL Scripts for Brain MRI Segmentation

## 1) Brain Extraction

Use `Bias field & neck cleanup` mode in the BET GUI and enter center coordinates for better brain extraction results.

```cmd
bet &
```

## 2) Linear Registration to MNI-152 Standard Space

```cmd
flirt -in tFUS33_brain.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -out tFUS33_T1_to_MNI_lin.nii.gz -omat tFUS33_T1_to_MNI_lin.mat
```

## 3) Non-Linear Registration

```cmd
fnirt --in=tFUS33_brain.nii.gz --aff=tFUS33_T1_to_MNI_lin.mat --ref=$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_1mm_brain_mask_dil.nii.gz --iout=tFUS33_T1_to_MNI_nonlin.nii.gz --cout=tFUS33_T1_to_MNI_warpcoef.nii.gz
```

## 4) Extract a Specific Region from the Juelich Atlas

We chose V5L in our case. You may choose your own ROI by checking the `$FSLDIR/data/atlases/Juelich.xml` for label index. 

*Note: +1 from the original index from the `.xml` file. For example, V5L: 88$\rightarrow$ 89 .*

```cmd
fslmaths $FSLDIR/data/atlases/Juelich/Juelich-maxprob-thr25-1mm.nii.gz -thr 89 -uthr 89 -bin V5_L_standard.nii.gz
```

Similarly, you can try V5R:

```cmd
fslmaths $FSLDIR/data/atlases/Juelich/Juelich-maxprob-thr25-1mm.nii.gz -thr 90 -uthr 90 -bin V5_R_standard.nii.gz
```

## 5) Invert the Transformation and Apply the Inverse Warp

```cmd
invwarp -w tFUS33_T1_to_MNI_warpcoef.nii.gz -o tFUS33_MNI_to_T1_warpcoef.nii.gz -r tFUS33_brain.nii.gz
```

### Generate binary masks (Case 1):

```cmd
applywarp --in=../V5_L_standard.nii.gz --ref=tFUS33_brain.nii.gz --out=tFUS33_V5_L_native.nii.gz --warp=tFUS33_MNI_to_T1_warpcoef.nii.gz
```

### Retain the T1 intensities (Case 2, for Kranion):

```cmd
applywarp --in=../V5_L_standard.nii.gz --ref=tFUS33_brain.nii.gz --warp=tFUS33_MNI_to_T1_warpcoef.nii.gz --out=tFUS33_V5_L_native.nii.gz --interp=nn

fslmaths tFUS33_brain.nii.gz -mul tFUS33_V5_L_native.nii.gz tFUS33_brain_V5_L_segmented.nii.gz
```
