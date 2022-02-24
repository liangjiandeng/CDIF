function showRGB8(Ori_HRMS,fusion,location)

th_MSrgb = image_quantile(Ori_HRMS(:,:,[5,3,2]), [0.01 0.995]);
I_fuse = image_stretch(fusion(:,:,[5,3,2]),th_MSrgb);
ent=rectangleonimage(I_fuse,location,0.5, 3,2, 2,4);
figure,imshow(ent,[])

end