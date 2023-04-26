clear;
clc;
message = 'choose the folder with images';
uiwait(msgbox(message));
root_dir = pwd;
dir_path = uigetdir();
cd(dir_path);
dir_contents = dir;
filenames = {dir_contents.name};
filenames = filenames(~ismember(filenames,{'.','..','.DS_Store','Thumbs.db'}));
loop = length(filenames)/5;
r = 0;
roi_set = cell(1,5);
for k = 0:(loop-1)
    indexs = [5*k+2,5*k+3,5*k+4,5*k+1];
    v=0;  %%while loop counter
    sl = 1;
    c=0;  %% filename counter
    nuclei = [];golgis = [];references = [];
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
    
    thetas = [0,0];distances = [0,0]; Areas = [0,0];Circularities = [0,0];
    Sh_indexes = [0,0];Perimeters = [0,0];
    nc = 0;
    while v == 0
        for i = 1:2
            %cell selection
            index = indexs(1);
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
            roi = [roi,{vert}]; %#ok<AGROW>
            fprintf(varfile,'Cell_set_%d\nx,y\n',sl);
            for j=1:size(vert,1)
                fprintf(varfile,'%f,%f\n',vert(j,:));
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
            Areas(i) = Area/(4.52*4.52);
            Area = Areas(i);
            Circularities(i) = s.Circularity;
            Perimeters(i) = s.Perimeter/4.52;
            Sh_indexes(i) = Perimeters(i)/sqrt(Area);
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
            golgi = [x_com_golgi,y_com_golgi];
            reference = [x_com_nucleus+50,y_com_nucleus];
            n2g = sqrt(sum(bsxfun(@minus, nucleus,golgi).^2,2))/4.52;
            distances(i) = n2g;
            x2 = x_com_golgi;
            y2 = y_com_golgi;
            x1 = x_com_nucleus;
            y1 = y_com_nucleus;
            if (x2-x1) < 0
                theta = pi - atan((y2-y1)/(x2-x1));
            else
                if (y2-y1) > 0
                    theta = 2*pi - atan((y2-y1)/(x2-x1));
                else
                    theta = -atan((y2-y1)/(x2-x1));
                end
            end
            thetas(i) = theta*180/pi;
            
            %         message=sprintf('distance in micrometer is %.2f \n angle in degrees is %.2f',distance,theta);
            %         uiwait(msgbox(message));
            %angle ends--------------------------------------------------------
            close all
           
            pause(2)
            
            nuclei = [nuclei;nucleus]; %#ok<AGROW>
            golgis = [golgis;golgi];    %#ok<AGROW>
            references = [references;reference]; %#ok<AGROW>
            nc = nc+1;
        end
        index =indexs(4);
        current_file = filenames(index);
        fullFileName = fullfile(dir_path,current_file);
        current_file = char(fullFileName);
        grayImage = imread(current_file);
        
        add = size(nuclei,1);
        nucleus = [nuclei(nc-1),nuclei(nc-1+add)];
        print = strcat(num2str(thetas(1)),' degree');
        text = insertText(grayImage,nucleus+1,print,'Fontsize',30);
        
        nucleus = [nuclei(nc),nuclei(nc+add)];
        print = strcat(num2str(thetas(2)),' degree');
        text = insertText(text,nucleus+1,print,'Fontsize',30);
        figure,imshow(text);
        hold on;
        axis on;
        for m = 0:1
            nucleus = [nuclei(m+nc-1),nuclei(m+nc-1+add)];
            golgi = [golgis(m+nc-1),golgis(m+nc-1+add)];
            reference = [references(m+nc-1);references(m+nc-1+add)];
            plot(nucleus(1),nucleus(2),'r.');
            plot(golgi(1),golgi(2),'y.');
            drawArrow(nucleus,golgi,'g'); %draws arrow from nucleus to golgi
            drawArrow(nucleus,reference,'g');
        end
        
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
            fprintf(fid,'sl,angle difference ,angle (cell_1), nucleus_golgi_distance(cell_1),angle (cell_2),nucleus_golgi_distance(cell_2),distance_golgis,distance_nuclei\n');
        end
        fullFileName = fullfile(folder,current_file);
        fullFileName = strcat(fullFileName,num2str(sl));
        %imwrite(gcf,fullFileName);
        saveas(gcf,fullFileName,'tif')
        close all
        
        difference = abs(thetas(1)-thetas(2));
        nucleus1 = [nuclei(nc-1),nuclei(nc+1)];
        nucleus2 = [nuclei(nc),nuclei(nc+2)];
        golgi1 = [golgis(nc-1),golgis(nc+1)];
        golgi2 = [golgis(nc),golgis(nc+2)];
        distance_golgis = sqrt(sum(bsxfun(@minus, golgi1,golgi2).^2,2))/4.52;
        distance_nuclei = sqrt(sum(bsxfun(@minus, nucleus1,nucleus2).^2,2))/4.52;
        fprintf(fid,'%d,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f\n',sl,difference,thetas(1),distances(1),thetas(2),distances(2),distance_golgis,distance_nuclei);
        
        if c==0
            message = sprintf('open an folder to save cell preperties data');
            uiwait(msgbox(message));
            dir = uigetdir();
        end
        name = sprintf('properties%d.csv',k);
        file = fullfile(dir,name);
        fid3 = fopen(file,'w');
        fprintf(fid3,'sl,Area,Perimeter,Circularity,Shape Index\n');
        
        
        if c==0
            message = sprintf('Now choose to use same file or move to next');
            uiwait(msgbox(message));
        end
        for n = 1:size(Areas,2)
            fprintf(fid3,'%d,%f,%f,%f,%f\n',n,Areas(n),Perimeters(n),Circularities(n),Sh_indexes(n));
        end
        v = input(' Enter 1 to move to next file or else enter 0 ');
        sl = sl+1;
        c = c+1;
    end
    roi_set(k+1) = {roi};
    fclose(fid);
    fclose(fid3);
    fclose(varfile);
    message = sprintf('open an folder to save nuclei,golgi and cell data');
    uiwait(msgbox(message));
    dir = uigetdir();

    name = sprintf('centers%d.csv',k);
    file = fullfile(dir,name);
    fid2 = fopen(file,'w');
    fprintf(fid2,'sl,nucleus_x,nucleus_y,golgi_x,golgi_y\n');
    x1 = nuclei(:,1);
    y1 = nuclei(:,2);
    x2 = golgis(:,1);
    y2 = golgis(:,2);
    for n = 1:size(nuclei,1)
        fprintf(fid2,'%d,%f,%f,%f,%f\n',n,x1(n),y1(n),x2(n),y2(n));
    end
    fclose(fid2);
end
cd(root_dir);