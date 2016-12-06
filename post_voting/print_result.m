dataroot    = getDataRoot;
res_img_dir = fullfile(dataroot,'post_result/201301270925');

result_file = fullfile(dataroot,'result/test1.mat');

load(result_file);

if(~exist(res_img_dir,'dir'))
    mkdir(res_img_dir);
end

for i=1:length(detData)
    fprintf(1, '.');
    img_file    = fullfile(dataroot,['images/', detData(i).imgname, '.png']);
    ax1 = subplot('Position', [0,0,1,1]);
    display_hypo_rect(imread(img_file), [], detData(i).hypo_list, detData(i).score_list, detData(i).bbox_list, ax1);
    axis image;
    print(gcf, '-dpng', fullfile(res_img_dir, [detData(i).imgname, '.png']));
    close(gcf);
end
fprintf(1,'\nDone\n');