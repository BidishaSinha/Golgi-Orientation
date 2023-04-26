clear;
clc;
message = 'choose the folder with images';
uiwait(msgbox(message));
dir_path = uigetdir(); 
dir_contents = dir(dir_path);
filenames = {dir_contents.name};
filenames = filenames(~ismember(filenames,{'.','..','.DS_Store'}));
loop = length(filenames)/5;
r = 0;
roi_set = cell(1,5);
for k = 0:4
    indexs = (5*k+1:5*(k+1));
    v=0;  %%while loop counter
    sl = 1;
    c=0;  %% filename counter
    nuclei = [];golgis = [];wound_points = [];
    %wound selection
    index = indexs(1);
    current_file = filenames(index);
    message = sprintf(' Draw boundary of the wound ');
    uiwait(msgbox(message));
    fullFileName = fullfile(dir_path,current_file);
    current_file = char(fullFileName);
    disp('The files name is : ' )
    disp(current_file)
    clf;
    Image = imread(current_file);    
    imshow(Image,[]);
    w = drawpolyline;
    pos = w.Position;
    close;
    N = size(Image,2);
    wound = interppolygon(pos,N);
    if r == 0
        message = sprintf(' select an folder to save your ROI data and also choose a namefor the file');
        uiwait(msgbox(message));
        dir_name = uigetdir();
        s1 = input(' Enter the name of ROI file');
    end
    s2 = sprintf('_Roi_%d.csv',k);
    roifile = strcat(s1,s2);
    roifile = fullfile(dir_name,roifile);
    varfile = fopen(roifile,'w');
    roi = {};
    r = r+1;
    while v == 0
        %cell selection
        index = indexs(5);
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        disp('The files name is : ' )
        disp(current_file)
        clf;
        Image = imread(current_file);    
        imshow(Image,[]);
        h = drawpolygon;
%         vert = h.Vertices;
        vert = h.Position;
        roi = [roi,{vert}];
        fprintf(varfile,'Cell_%d\nx,y\n',sl);
        for i=1:size(vert,1)
            fprintf(varfile,'%f,%f\n',vert(i,:));
        end
        original_mask = h.createMask();
        mask = uint16(original_mask);
%         imshow(original_mask);
        index = indexs(1);
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        disp('The files name is : ' )
        disp(current_file)
        Image = imread(current_file);
        matrix = Image.*mask;
        X = imbinarize(matrix,graythresh(matrix));
        s = regionprops(X,'Area','Circularity','Perimeter');
        Area = s.Area;
        Area = Area/(6.022*6.022);
        Circularity = s.Circularity;
        Perimeter = s.Perimeter/6.022;
        sh.index = Perimeter/sqrt(Area);
         % golgi----------------------------------------------------------
         % selectoin------------------------------------------------------=
        index = indexs(2);       
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        disp('The files name is : ' )
        disp(current_file)
        grayImage = imread(current_file);
        imshow(grayImage,[]);
        hold on
%         mask = original_mask;
%         mask = uint16(mask);
        matrix = grayImage.*mask;
        r =  max(matrix(:))*0.4 ;
        matrix(matrix <r) = 0;
        original_picture_matrix = matrix;
        A = matrix;
        A(A > 0 ) = 1;
        binarised_picture_matrix = A;
        A = double(original_picture_matrix);
        B = double(binarised_picture_matrix);
        [ii,jj] = ndgrid(1:size(A,1),1:size(A,2));
        y_com_golgi = sum(sum(ii.*(A.*B)))/sum(sum(A.*B));
        x_com_golgi = sum(sum(jj.*(A.*B)))/sum(sum(A.*B));
        disp(x_com_golgi);disp(y_com_golgi);
        plot(x_com_golgi,y_com_golgi,'r*');
        pause(2);
        clf;
        %golgi selection ends----------------------------------------------
        %Nucleus Selection-------------------------------------------------
        index = indexs(3);
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        disp('The files name is : ')
        disp(current_file)
        grayImage = imread(current_file);    
        imshow(grayImage,[]);
        hold on
%         mask = h.createMask();
        A = grayImage;
        B = mask;
        A = A.*B;
        r =  max(A(:))*0.4 ;
        A(A <r) = 0;
        B = double(A.*B);
        [ii,jj] = ndgrid(1:size(A,1),1:size(A,2));
        y_com_nucleus = sum(ii.*B)/sum(B);
        x_com_nucleus = sum(jj.*B)/sum(B);
        disp(x_com_nucleus);disp(y_com_nucleus);
        imshow(grayImage,[]);
        hold on;
        plot(x_com_nucleus,y_com_nucleus,'r*');
        pause(2)
        %Nuleus selection ends---------------------------------------------
        %Angle Calculation-------------------------------------------------
        nucleus = [x_com_nucleus, y_com_nucleus];
        distances = sqrt(sum(bsxfun(@minus, nucleus,wound).^2,2));
        distance = min(distances);
        closest = wound(find(distances==min(distances)),:);  %#ok<*FNDSB>
        golgi = [x_com_golgi,y_com_golgi];
        wound_point = closest;
        n2w = sum(bsxfun(@minus, nucleus,wound_point).^2,2);
        n2g = sum(bsxfun(@minus, nucleus,golgi).^2,2);
        w2g = sum(bsxfun(@minus, golgi,wound_point).^2,2);
        theta = acos((n2w+n2g-w2g)/(2*sqrt(n2w*n2g)));
        theta = theta*180/pi;
        
%         message=sprintf('distance in micrometer is %.2f \n angle in degrees is %.2f',distance,theta);
%         uiwait(msgbox(message));
        %angle ends--------------------------------------------------------
        close all
        index =indexs(4);
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        grayImage = imread(current_file);
        print = strcat(num2str(theta),' degree');
        text = insertText(grayImage,nucleus+1,print,'Fontsize',30);
        figure,imshow(text);
        hold on;
        axis on;
        plot(wound(:,1),wound(:,2),'r');
        plot(nucleus(1),nucleus(2),'r.'); 
        plot(golgi(1),golgi(2),'y.');
        plot(wound_point(1),wound_point(2),'b.');
        drawArrow(nucleus,golgi,'g'); %draws arrow from nucleus to golgi
        drawArrow(nucleus,wound_point,'g'); %draws arrow from nucleus to wound
%         plot(vert(:,1),vert(:,2));
        pause(2)
        current_file = char(filenames(index));
        l = length(current_file);
        current_file(l-3:l) = '';
        if c == 0
            message = sprintf(' select an folder to save your image ');
            uiwait(msgbox(message));
            folder = uigetdir();
            results = strcat(current_file,'.csv');
            fileName = fullfile(folder, results);
            fid = fopen(fileName, 'wt');
            fprintf(fid,'sl,angle(degrees),distance(cell to wound), Area, Circularity, Perimeter, Shape Index\n');
        end
        fprintf(fid,'%d,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%f\n',sl,theta,distance,Area,Circularity,Perimeter,sh.index);
        fullFileName = fullfile(folder,current_file);
        fullFileName = strcat(fullFileName,num2str(sl));
        %imwrite(gcf,fullFileName);
        saveas(gcf,fullFileName,'tif')
        close all
        if c==0
            message = sprintf('Now choose to use same file or move to next');
            uiwait(msgbox(message));
        end
        v = input(' Enter 1 to move to next file or else enter 0 ');
        sl = sl+1;
        c = c+1;
        nuclei = [nuclei;nucleus]; %#ok<AGROW>
        golgis = [golgis;golgi];    %#ok<AGROW>
        wound_points = [wound_points;wound_point];  %#ok<AGROW>
        
    end
    roi_set(k+1) = {roi};
    fclose(fid);
    fclose(varfile);
    message = sprintf('open an folder to save nuclei,golgi and wound data');
    uiwait(msgbox(message));
    dir = uigetdir();

    name = sprintf('centers%d.csv',k);
    file = fullfile(dir,name);
    fid2 = fopen(file,'w');
    fprintf(fid2,'sl,nucleus_x,nucleus_y,golgi_x,golgi_y,wound_x,wound_y\n');
    x3 = wound_points(:,1);
    y3 = wound_points(:,2);
    x1 = nuclei(:,1);
    y1 = nuclei(:,2);
    x2 = golgis(:,1);
    y2 = golgis(:,2);
    for n = 1:size(nuclei,1)
        fprintf(fid2,'%d,%f,%f,%f,%f,%f,%f\n',n,x1(n),y1(n),x2(n),y2(n),x3(n),y3(n));
    end
    fclose(fid2);
end
