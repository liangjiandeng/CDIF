function showRGB4(Ori_HRMS,fusion,location)

th_MSrgb = image_quantile(Ori_HRMS(:,:,[3,2,1]), [0.01 0.99]);
I_fuse = image_stretch(fusion(:,:,[3,2,1]),th_MSrgb);
ent=rectangleonimage(I_fuse,location, 0.5, 3, 1, 2, 1);
figure,imshow(ent,[])

end