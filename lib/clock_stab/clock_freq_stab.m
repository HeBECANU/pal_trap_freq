% process video/ pictures of frequ counter to get out the freq stability

%vid_file='F:\WIN_20220117_17_46_12_Pro.mp4'
vid_file='F:\WIN_20220117_17_46_42_Pro.mp4';
v=VideoReader(vid_file)

%%

max_frames=floor(v.FrameRate*v.Duration);
time_between_poll=2;
frames_between_poll=ceil(v.FrameRate*time_between_poll);
time_between_poll=frames_between_poll/v.FrameRate;


%%
test_img=read(v,100);
rotation_angle=2;
test_img=imrotate(test_img,rotation_angle, 'bicubic');
figure(2)
imshow(test_img)
%roi = drawrectangle('Color','g');
% roi_pos=round(roi.Position);
roi_pos=[567   528   655    87];
im_roi=test_img(roi_pos(2):(roi_pos(2)+roi_pos(4)),roi_pos(1):(roi_pos(1)+roi_pos(3)),:);
imshow(im_roi)


%%
%poll_frames_idx=1:frames_between_poll:(frames_between_poll*100);
poll_frames_idx=1:frames_between_poll:max_frames;
iimax=numel(poll_frames_idx);
poll_frames_img=cell(iimax,1);
fprintf('getting frames')
for ii=1:iimax
    fprintf('%03uof%03u\n',ii,iimax)
    im_tmp=read(v,poll_frames_idx(ii));
    im_tmp=imrotate(im_tmp,rotation_angle, 'bicubic');
    im_tmp=im_tmp(roi_pos(2):(roi_pos(2)+roi_pos(4)),roi_pos(1):(roi_pos(1)+roi_pos(3)),:);
    poll_frames_img{ii}=im_tmp;
    %imwrite(im_tmp,sprintf('F:\\subregions\\%u.jpg',ii))
end

%% run ocr
trainedLanguage = 'E:\Dropbox\UNI\programs\Trap_freq_methods\lib\freq_count_ocr\sevensegment\tessdata\sevensegment.traineddata';

frames_text=cell(iimax,1);
frames_ocr_obj=cell(iimax,1);
iimax=numel(poll_frames_img);
fprintf('running ocr')
for ii=1:iimax
    fprintf('%03uof%03u\n',ii,iimax)
    im_ocr_tmp=poll_frames_img{ii};
    %im_roi_binarized=im_roi(:,:,2)>190;
    ocrResults = ocr(im_ocr_tmp,'CharacterSet','0123456789.','Language', trainedLanguage, ...
            'TextLayout', 'Block');
    frames_ocr_obj{ii}=ocrResults;
    frames_text{ii}=strtrim(ocrResults.Text);
end

%% process into a vector of values
freq_meas=[];
freq_meas.times=poll_frames_idx/v.FrameRate;
freq_meas.freqs=nan(iimax,1)
iimax=numel(freq_meas.times);
fprintf('converting to vector')
for ii=1:iimax
    str_tmp=frames_text{ii};
    str_tmp=strrep(str_tmp,' ','')
    if str_tmp(end)~='3'
        warning('not ending with 3')
        pause
    end
    if str_tmp(3)~='.'
        warning('does not have a . inserting')
        str_tmp=insertAfter(str_tmp,2,'.');
    end
    %imshow(poll_frames_img{ii})
    if numel(str_tmp)~=11
        freq_meas.freqs(ii)=nan;
    else
        value=str2num(strtrim(str_tmp(1:end-1)))*1e3;
        freq_meas.freqs(ii)=value;
    end
    
end

%%
window_size_sec=60;
window_size_idx=round(window_size_sec/mean(diff(freq_meas.times)));
freq_outlier=isoutlier(freq_meas.freqs,'movmedian',window_size_idx,'ThresholdFactor',4);
freq_meas.freqs(freq_outlier)=nan;

freq_meas.freqs=fillmissing(freq_meas.freqs,'pchip');

%%
freq_mean=mean(freq_meas.freqs);
freq_std=std(freq_meas.freqs);
freq_frac_dev=(freq_meas.freqs-freq_mean)/freq_mean
plot(freq_meas.times,freq_frac_dev)


%% 
%figure(1)
%imshow(imrotate(im_roi,2))

%%
% im_ocr_test=poll_frames_img{4};
% trainedLanguage = 'E:\Dropbox\UNI\programs\Trap_freq_methods\lib\freq_count_ocr\sevensegment\tessdata\sevensegment.traineddata';
% 
% %im_roi_binarized=im_roi(:,:,2)>190;
% ocrResults = ocr(im_ocr_test,'CharacterSet','0123456789.','Language', trainedLanguage, ...
%         'TextLayout', 'Block')
% regularExpr='\d';
% bboxes = locateText(ocrResults, regularExpr,'UseRegexp',true);
% digits = regexp(ocrResults.Text,regularExpr,'match');
% Iocr = insertObjectAnnotation(im_ocr_test,'rectangle',bboxes,digits);
% figure(1);
% imshow(Iocr);



% %% training data
% ocr_train_img=zeros(200,2000,3);
% train_font=repmat('1 2 3 4 5 6 7 8 9 0 . ',[1,2])
% ocr_train_img = insertText(ocr_train_img,[0 0],train_font,...
%     'Font','Seven Segment','FontSize',100,'BoxOpacity',0,'TextColor','w');
% H = fspecial('disk',5);
% ocr_train_img = imfilter(ocr_train_img,H,'replicate'); 
% ocr_train_img=sum(ocr_train_img,3)>0.8;
% imshow(ocr_train_img)
