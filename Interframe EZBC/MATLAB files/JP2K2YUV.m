clear all;
clc;

num_frames = 144;
GOPsz = 32;

GOPnum = num_frames / GOPsz;
GOPres = num_frames - GOPsz * floor(GOPnum);
GOPnum = floor(GOPnum);

width = 352;
height = 288;
y_size = width*height;
c_size = y_size / 4;
f_size = y_size + 2*c_size;

Y = zeros(y_size*num_frames,1);
U = zeros(c_size*num_frames,1);
V = zeros(c_size*num_frames,1);

outfid = fopen('BUS_352x288_30_dec_jp2k_01.yuv','wb');

% %expand JP2K data and write back into subbands
% for i = 0:(GOPnum-1)
% 
%      jpx_name = ['tmp_mob' num2str(i) '.j2c'];
% %     system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl ' ' -raw_components ']);
%      system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,y64.rawl,y65.rawl,y66.rawl,y67.rawl,y68.rawl,y69.rawl,y70.rawl,y71.rawl,y72.rawl,y73.rawl,y74.rawl,y75.rawl,y76.rawl,y77.rawl,y78.rawl,y79.rawl,y80.rawl,y81.rawl,y82.rawl,y83.rawl,y84.rawl,y85.rawl,y86.rawl,y87.rawl,y88.rawl,y89.rawl,y90.rawl,y91.rawl,y92.rawl,y93.rawl,y94.rawl,y95.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,u32.rawl,u33.rawl,u34.rawl,u35.rawl,u36.rawl,u37.rawl,u38.rawl,u39.rawl,u40.rawl,u41.rawl,u42.rawl,u43.rawl,u44.rawl,u45.rawl,u46.rawl,u47.rawl,u48.rawl,u49.rawl,u50.rawl,u51.rawl,u52.rawl,u53.rawl,u54.rawl,u55.rawl,u56.rawl,u57.rawl,u58.rawl,u59.rawl,u60.rawl,u61.rawl,u62.rawl,u63.rawl,u64.rawl,u65.rawl,u66.rawl,u67.rawl,u68.rawl,u69.rawl,u70.rawl,u71.rawl,u72.rawl,u73.rawl,u74.rawl,u75.rawl,u76.rawl,u77.rawl,u78.rawl,u79.rawl,u80.rawl,u81.rawl,u82.rawl,u83.rawl,u84.rawl,u85.rawl,u86.rawl,u87.rawl,u88.rawl,u89.rawl,u90.rawl,u91.rawl,u92.rawl,u93.rawl,u94.rawl,u95.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl,v32.rawl,v33.rawl,v34.rawl,v35.rawl,v36.rawl,v37.rawl,v38.rawl,v39.rawl,v40.rawl,v41.rawl,v42.rawl,v43.rawl,v44.rawl,v45.rawl,v46.rawl,v47.rawl,v48.rawl,v49.rawl,v50.rawl,v51.rawl,v52.rawl,v53.rawl,v54.rawl,v55.rawl,v56.rawl,v57.rawl,v58.rawl,v59.rawl,v60.rawl,v61.rawl,v62.rawl,v63.rawl,v64.rawl,v65.rawl,v66.rawl,v67.rawl,v68.rawl,v69.rawl,v70.rawl,v71.rawl,v72.rawl,v73.rawl,v74.rawl,v75.rawl,v76.rawl,v77.rawl,v78.rawl,v79.rawl,v80.rawl,v81.rawl,v82.rawl,v83.rawl,v84.rawl,v85.rawl,v86.rawl,v87.rawl,v88.rawl,v89.rawl,v90.rawl,v91.rawl,v92.rawl,v93.rawl,v94.rawl,v95.rawl ' ' -raw_components ']);
% 
% %    system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,u32.rawl,u33.rawl,u34.rawl,u35.rawl,u36.rawl,u37.rawl,u38.rawl,u39.rawl,u40.rawl,u41.rawl,u42.rawl,u43.rawl,u44.rawl,u45.rawl,u46.rawl,u47.rawl,u48.rawl,u49.rawl,u50.rawl,u51.rawl,u52.rawl,u53.rawl,u54.rawl,u55.rawl,u56.rawl,u57.rawl,u58.rawl,u59.rawl,u60.rawl,u61.rawl,u62.rawl,u63.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl,v32.rawl,v33.rawl,v34.rawl,v35.rawl,v36.rawl,v37.rawl,v38.rawl,v39.rawl,v40.rawl,v41.rawl,v42.rawl,v43.rawl,v44.rawl,v45.rawl,v46.rawl,v47.rawl,v48.rawl,v49.rawl,v50.rawl,v51.rawl,v52.rawl,v53.rawl,v54.rawl,v55.rawl,v56.rawl,v57.rawl,v58.rawl,v59.rawl,v60.rawl,v61.rawl,v62.rawl,v63.rawl ' ' -raw_components ']);
% %     system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,y64.rawl,y65.rawl,y66.rawl,y67.rawl,y68.rawl,y69.rawl,y70.rawl,y71.rawl,y72.rawl,y73.rawl,y74.rawl,y75.rawl,y76.rawl,y77.rawl,y78.rawl,y79.rawl,y80.rawl,y81.rawl,y82.rawl,y83.rawl,y84.rawl,y85.rawl,y86.rawl,y87.rawl,y88.rawl,y89.rawl,y90.rawl,y91.rawl,y92.rawl,y93.rawl,y94.rawl,y95.rawl ' ' -raw_components ']);
% 
%     for j = 0:(GOPsz-1)
%   
%         fid_name1 = ['y' num2str( j ) '.rawl'];
%         fid_name2 = ['u' num2str( j ) '.rawl']; 
%         fid_name3 = ['v' num2str( j ) '.rawl']; 
% 
%         fp = fopen(fid_name1,'rb');
%         Y = fread(fp,y_size,'int16');
%         fclose(fp);
% 
%         fp = fopen(fid_name2,'rb');
%         U = fread(fp,c_size,'int16');
%         fclose(fp);
% 
%         fp = fopen(fid_name3,'rb');
%         V = fread(fp,c_size,'int16');
%         fclose(fp);
% 
% %         fp = fopen(fid_name1,'rb');
% %         Y = fread(fp,y_size,'int16');
% %         U = fread(fp,c_size,'int16');
% %         V = fread(fp,c_size,'int16');
% %         fclose(fp);
% 
%         Y = (Y/32);
%         U = (U/32);
%         V = (V/32);
%         
% %         if j == 0 || j == 32
% %             U = U + 128;
% %             V = V + 128;
% %         end
% 
%         fwrite(outfid,Y,'float32');
%         fwrite(outfid,U,'float32');
%         fwrite(outfid,V,'float32');
%         
% %         if j == 0
% %             Yd = reshape(Y,width,height);
% %             Ud = reshape(U,width/2,height/2);
% %             Vd = reshape(V,width/2,height/2);      
% %         end
%         
%     end
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOPsz = 144;
% GOPnum = 1;
% GOPres = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%expand JP2K data and write back into subbands
for i = 0:(GOPnum-1)

     jpx_name = ['tmp_bq' num2str(i) '.j2c'];
%     system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl ' ' -raw_components ']);
%     system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,y64.rawl,y65.rawl,y66.rawl,y67.rawl,y68.rawl,y69.rawl,y70.rawl,y71.rawl,y72.rawl,y73.rawl,y74.rawl,y75.rawl,y76.rawl,y77.rawl,y78.rawl,y79.rawl,y80.rawl,y81.rawl,y82.rawl,y83.rawl,y84.rawl,y85.rawl,y86.rawl,y87.rawl,y88.rawl,y89.rawl,y90.rawl,y91.rawl,y92.rawl,y93.rawl,y94.rawl,y95.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,u32.rawl,u33.rawl,u34.rawl,u35.rawl,u36.rawl,u37.rawl,u38.rawl,u39.rawl,u40.rawl,u41.rawl,u42.rawl,u43.rawl,u44.rawl,u45.rawl,u46.rawl,u47.rawl,u48.rawl,u49.rawl,u50.rawl,u51.rawl,u52.rawl,u53.rawl,u54.rawl,u55.rawl,u56.rawl,u57.rawl,u58.rawl,u59.rawl,u60.rawl,u61.rawl,u62.rawl,u63.rawl,u64.rawl,u65.rawl,u66.rawl,u67.rawl,u68.rawl,u69.rawl,u70.rawl,u71.rawl,u72.rawl,u73.rawl,u74.rawl,u75.rawl,u76.rawl,u77.rawl,u78.rawl,u79.rawl,u80.rawl,u81.rawl,u82.rawl,u83.rawl,u84.rawl,u85.rawl,u86.rawl,u87.rawl,u88.rawl,u89.rawl,u90.rawl,u91.rawl,u92.rawl,u93.rawl,u94.rawl,u95.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl,v32.rawl,v33.rawl,v34.rawl,v35.rawl,v36.rawl,v37.rawl,v38.rawl,v39.rawl,v40.rawl,v41.rawl,v42.rawl,v43.rawl,v44.rawl,v45.rawl,v46.rawl,v47.rawl,v48.rawl,v49.rawl,v50.rawl,v51.rawl,v52.rawl,v53.rawl,v54.rawl,v55.rawl,v56.rawl,v57.rawl,v58.rawl,v59.rawl,v60.rawl,v61.rawl,v62.rawl,v63.rawl,v64.rawl,v65.rawl,v66.rawl,v67.rawl,v68.rawl,v69.rawl,v70.rawl,v71.rawl,v72.rawl,v73.rawl,v74.rawl,v75.rawl,v76.rawl,v77.rawl,v78.rawl,v79.rawl,v80.rawl,v81.rawl,v82.rawl,v83.rawl,v84.rawl,v85.rawl,v86.rawl,v87.rawl,v88.rawl,v89.rawl,v90.rawl,v91.rawl,v92.rawl,v93.rawl,v94.rawl,v95.rawl ' ' -raw_components ']);

%    system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,u32.rawl,u33.rawl,u34.rawl,u35.rawl,u36.rawl,u37.rawl,u38.rawl,u39.rawl,u40.rawl,u41.rawl,u42.rawl,u43.rawl,u44.rawl,u45.rawl,u46.rawl,u47.rawl,u48.rawl,u49.rawl,u50.rawl,u51.rawl,u52.rawl,u53.rawl,u54.rawl,u55.rawl,u56.rawl,u57.rawl,u58.rawl,u59.rawl,u60.rawl,u61.rawl,u62.rawl,u63.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl,v32.rawl,v33.rawl,v34.rawl,v35.rawl,v36.rawl,v37.rawl,v38.rawl,v39.rawl,v40.rawl,v41.rawl,v42.rawl,v43.rawl,v44.rawl,v45.rawl,v46.rawl,v47.rawl,v48.rawl,v49.rawl,v50.rawl,v51.rawl,v52.rawl,v53.rawl,v54.rawl,v55.rawl,v56.rawl,v57.rawl,v58.rawl,v59.rawl,v60.rawl,v61.rawl,v62.rawl,v63.rawl ' ' -raw_components ']);
%     system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,y64.rawl,y65.rawl,y66.rawl,y67.rawl,y68.rawl,y69.rawl,y70.rawl,y71.rawl,y72.rawl,y73.rawl,y74.rawl,y75.rawl,y76.rawl,y77.rawl,y78.rawl,y79.rawl,y80.rawl,y81.rawl,y82.rawl,y83.rawl,y84.rawl,y85.rawl,y86.rawl,y87.rawl,y88.rawl,y89.rawl,y90.rawl,y91.rawl,y92.rawl,y93.rawl,y94.rawl,y95.rawl ' ' -raw_components ']);

%    system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,y32.rawl,y33.rawl,y34.rawl,y35.rawl,y36.rawl,y37.rawl,y38.rawl,y39.rawl,y40.rawl,y41.rawl,y42.rawl,y43.rawl,y44.rawl,y45.rawl,y46.rawl,y47.rawl,y48.rawl,y49.rawl,y50.rawl,y51.rawl,y52.rawl,y53.rawl,y54.rawl,y55.rawl,y56.rawl,y57.rawl,y58.rawl,y59.rawl,y60.rawl,y61.rawl,y62.rawl,y63.rawl,y64.rawl,y65.rawl,y66.rawl,y67.rawl,y68.rawl,y69.rawl,y70.rawl,y71.rawl,y72.rawl,y73.rawl,y74.rawl,y75.rawl,y76.rawl,y77.rawl,y78.rawl,y79.rawl,y80.rawl,y81.rawl,y82.rawl,y83.rawl,y84.rawl,y85.rawl,y86.rawl,y87.rawl,y88.rawl,y89.rawl,y90.rawl,y91.rawl,y92.rawl,y93.rawl,y94.rawl,y95.rawl,y96.rawl,y97.rawl,y98.rawl,y99.rawl,y100.rawl,y101.rawl,y102.rawl,y103.rawl,y104.rawl,y105.rawl,y106.rawl,y107.rawl,y108.rawl,y109.rawl,y110.rawl,y111.rawl,y112.rawl,y113.rawl,y114.rawl,y115.rawl,y116.rawl,y117.rawl,y118.rawl,y119.rawl,y120.rawl,y121.rawl,y122.rawl,y123.rawl,y124.rawl,y125.rawl,y126.rawl,y127.rawl,y128.rawl,y129.rawl,y130.rawl,y131.rawl,y132.rawl,y133.rawl,y134.rawl,y135.rawl,y136.rawl,y137.rawl,y138.rawl,y139.rawl,y140.rawl,y141.rawl,y142.rawl,y143.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,u32.rawl,u33.rawl,u34.rawl,u35.rawl,u36.rawl,u37.rawl,u38.rawl,u39.rawl,u40.rawl,u41.rawl,u42.rawl,u43.rawl,u44.rawl,u45.rawl,u46.rawl,u47.rawl,u48.rawl,u49.rawl,u50.rawl,u51.rawl,u52.rawl,u53.rawl,u54.rawl,u55.rawl,u56.rawl,u57.rawl,u58.rawl,u59.rawl,u60.rawl,u61.rawl,u62.rawl,u63.rawl,u64.rawl,u65.rawl,u66.rawl,u67.rawl,u68.rawl,u69.rawl,u70.rawl,u71.rawl,u72.rawl,u73.rawl,u74.rawl,u75.rawl,u76.rawl,u77.rawl,u78.rawl,u79.rawl,u80.rawl,u81.rawl,u82.rawl,u83.rawl,u84.rawl,u85.rawl,u86.rawl,u87.rawl,u88.rawl,u89.rawl,u90.rawl,u91.rawl,u92.rawl,u93.rawl,u94.rawl,u95.rawl,u96.rawl,u97.rawl,u98.rawl,u99.rawl,u100.rawl,u101.rawl,u102.rawl,u103.rawl,u104.rawl,u105.rawl,u106.rawl,u107.rawl,u108.rawl,u109.rawl,u110.rawl,u111.rawl,u112.rawl,u113.rawl,u114.rawl,u115.rawl,u116.rawl,u117.rawl,u118.rawl,u119.rawl,u120.rawl,u121.rawl,u122.rawl,u123.rawl,u124.rawl,u125.rawl,u126.rawl,u127.rawl,u128.rawl,u129.rawl,u130.rawl,u131.rawl,u132.rawl,u133.rawl,u134.rawl,u135.rawl,u136.rawl,u137.rawl,u138.rawl,u139.rawl,u140.rawl,u141.rawl,u142.rawl,u143.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl,v32.rawl,v33.rawl,v34.rawl,v35.rawl,v36.rawl,v37.rawl,v38.rawl,v39.rawl,v40.rawl,v41.rawl,v42.rawl,v43.rawl,v44.rawl,v45.rawl,v46.rawl,v47.rawl,v48.rawl,v49.rawl,v50.rawl,v51.rawl,v52.rawl,v53.rawl,v54.rawl,v55.rawl,v56.rawl,v57.rawl,v58.rawl,v59.rawl,v60.rawl,v61.rawl,v62.rawl,v63.rawl,v64.rawl,v65.rawl,v66.rawl,v67.rawl,v68.rawl,v69.rawl,v70.rawl,v71.rawl,v72.rawl,v73.rawl,v74.rawl,v75.rawl,v76.rawl,v77.rawl,v78.rawl,v79.rawl,v80.rawl,v81.rawl,v82.rawl,v83.rawl,v84.rawl,v85.rawl,v86.rawl,v87.rawl,v88.rawl,v89.rawl,v90.rawl,v91.rawl,v92.rawl,v93.rawl,v94.rawl,v95.rawl,v96.rawl,v97.rawl,v98.rawl,v99.rawl,v100.rawl,v101.rawl,v102.rawl,v103.rawl,v104.rawl,v105.rawl,v106.rawl,v107.rawl,v108.rawl,v109.rawl,v110.rawl,v111.rawl,v112.rawl,v113.rawl,v114.rawl,v115.rawl,v116.rawl,v117.rawl,v118.rawl,v119.rawl,v120.rawl,v121.rawl,v122.rawl,v123.rawl,v124.rawl,v125.rawl,v126.rawl,v127.rawl,v128.rawl,v129.rawl,v130.rawl,v131.rawl,v132.rawl,v133.rawl,v134.rawl,v135.rawl,v136.rawl,v137.rawl,v138.rawl,v139.rawl,v140.rawl,v141.rawl,v142.rawl,v143.rawl ' ' -raw_components ']);


    system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,y16.rawl,y17.rawl,y18.rawl,y19.rawl,y20.rawl,y21.rawl,y22.rawl,y23.rawl,y24.rawl,y25.rawl,y26.rawl,y27.rawl,y28.rawl,y29.rawl,y30.rawl,y31.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,u16.rawl,u17.rawl,u18.rawl,u19.rawl,u20.rawl,u21.rawl,u22.rawl,u23.rawl,u24.rawl,u25.rawl,u26.rawl,u27.rawl,u28.rawl,u29.rawl,u30.rawl,u31.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl,v16.rawl,v17.rawl,v18.rawl,v19.rawl,v20.rawl,v21.rawl,v22.rawl,v23.rawl,v24.rawl,v25.rawl,v26.rawl,v27.rawl,v28.rawl,v29.rawl,v30.rawl,v31.rawl ' ' -raw_components ']);
    
    for k = 0:(GOPsz-1)
            fid_name1 = ['y' num2str( k ) '.rawl'];
            fid_name2 = ['u' num2str( k ) '.rawl']; 
            fid_name3 = ['v' num2str( k ) '.rawl']; 

            fp = fopen(fid_name1,'rb');
            Y = fread(fp,y_size,'int16');
            fclose(fp);

            fp = fopen(fid_name2,'rb');
            U = fread(fp,c_size,'int16');
            fclose(fp);

            fp = fopen(fid_name3,'rb');
            V = fread(fp,c_size,'int16');
            fclose(fp);

            Y = (Y/16);
            U = (U/16);
            V = (V/16);

            fwrite(outfid,Y,'float32');
            fwrite(outfid,U,'float32');
            fwrite(outfid,V,'float32');
    end     
    
end


if GOPres ~= 0
    jpx_name = ('tmp_bq_res.j2c');
    system([ ' kdu_expand ' ' -i ' jpx_name ' -o ' ' y0.rawl,y1.rawl,y2.rawl,y3.rawl,y4.rawl,y5.rawl,y6.rawl,y7.rawl,y8.rawl,y9.rawl,y10.rawl,y11.rawl,y12.rawl,y13.rawl,y14.rawl,y15.rawl,u0.rawl,u1.rawl,u2.rawl,u3.rawl,u4.rawl,u5.rawl,u6.rawl,u7.rawl,u8.rawl,u9.rawl,u10.rawl,u11.rawl,u12.rawl,u13.rawl,u14.rawl,u15.rawl,v0.rawl,v1.rawl,v2.rawl,v3.rawl,v4.rawl,v5.rawl,v6.rawl,v7.rawl,v8.rawl,v9.rawl,v10.rawl,v11.rawl,v12.rawl,v13.rawl,v14.rawl,v15.rawl ' ' -raw_components ']);

    for k = 0:(GOPres-1)
            fid_name1 = ['y' num2str( k ) '.rawl'];
            fid_name2 = ['u' num2str( k ) '.rawl']; 
            fid_name3 = ['v' num2str( k ) '.rawl']; 

            fp = fopen(fid_name1,'rb');
            Y = fread(fp,y_size,'int16');
            fclose(fp);

            fp = fopen(fid_name2,'rb');
            U = fread(fp,c_size,'int16');
            fclose(fp);

            fp = fopen(fid_name3,'rb');
            V = fread(fp,c_size,'int16');
            fclose(fp);

            Y = (Y/16);
            U = (U/16);
            V = (V/16);

            fwrite(outfid,Y,'float32');
            fwrite(outfid,U,'float32');
            fwrite(outfid,V,'float32');
    end     

end

fclose(outfid);