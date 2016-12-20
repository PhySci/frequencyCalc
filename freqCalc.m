% Class for processing results of OOMMF simulations 
classdef freqCalc < hgsetget % subclass hgsetget
   
properties
    meshunit = 'm';
    meshtype = 'rectangular';
    xbase
    ybase
    zbase
    xnodes
    ynodes
    znodes
    xstepsize = 0.001;
    ystepsize = 0.001;
    zstepsize = 0.001;
    xmin = 0;
    ymin = 0;
    zmin = 0;
    xmax = 0.01;
    ymax = 0.01;
    zmax = 0.01;
    dim = 3;
    
    % ground magnetization
    Mx0;
    My0;
    Mz0;
    Ms;
    
    % dynamic fields in local coordinate system
    hInpLoc = [];
    hOutLoc = [];
    
    % demagnetizing factors
    Nxx = [];
    Nyy = [];
    Nzz = [];
    Nyx = [];
    Nxy = [];
    
    % container for array of magnetic susceptibility tensor
    suscept = containers.Map('KeyType','double','ValueType','any');
    
    % FMR frequency
    freq = [];
    
    % external magentic field (Oe)
    H = 530;
    
    % damping constant
    alpha = 0.001
    
    % conversion factor for magnetization (A/m -> G)
    convM = 4*pi/1e3;
    % conversion factor for magnetic field (T -> G)
    convH = 1e4;
    
    % giromagnetic ratio
    gamma = 1.76e7;
end

properties (Access = private)
    % name of files
    
    % magnetization files
    staticFile = 'static.stc';
    outFile = 'outplane.stc';
    inpFile = 'inplane.stc';
    
    % effective field files
    BeffFile = 'B_eff_static.ovf'; 
    BeffOutFile = 'B_eff_outplane.ovf'; 
    BeffInpFile = 'B_eff_inpplane.ovf'; 

    % demagnetizing field files
    BdemagInpFile = 'B_demag_inplane.ovf';
    BdemagOutFile = 'B_demag_outplane.ovf';
    BdemagStaticFile = 'B_demag_static.ovf';
   
    % exchange field files
    BexchInpFile = 'B_exch_inplane.ovf';
    BexchOutFile = 'B_exch_outplane.ovf';
    BexchStaticFile = 'B_exch_static.ovf';
    
    % deflection amplitude
    def = 1e-4;
        
    header = '';
    
    format = '';
end 

properties (Constant)
    doubleTestVal = 123456789012345.0;
    singleTestVal = 1234567.0;
end

methods
    
    % constructor
    function obj = freqCalc()
    end
    
    % calculate in-plane and out-of-plane deflections of magnetisation,
    % save deflected states
    function calcDeflection(obj)
        obj.readFile(obj.staticFile);
        disp('Calculate fields');
        M0 = obj.readFile(obj.staticFile);
        Ms = sqrt(M0(:,:,:,1).^2+M0(:,:,:,2).^2+M0(:,:,:,3).^2);
        
        % calculate direction of in-plane deflection
        norm = Ms./sqrt(M0(:,:,:,1).^2+M0(:,:,:,2).^2);
        nInp(:,:,:,1) = -norm.*M0(:,:,:,2);
        nInp(:,:,:,2) = norm.*M0(:,:,:,1);
        nInp(:,:,:,3) = zeros(obj.xnodes,obj.ynodes,obj.znodes);
        Minp = M0 + obj.def*nInp;
        
        % calculate direction of out-plane deflection
        norm = sqrt(M0(:,:,:,1).^2+M0(:,:,:,2).^2); 
        nOut(:,:,:,1) = (-M0(:,:,:,1).*M0(:,:,:,3))./norm;
        nOut(:,:,:,2) = (-M0(:,:,:,2).*M0(:,:,:,3))./norm;
        nOut(:,:,:,3) = norm;
        Mout = M0 + obj.def*nOut;      
        
        disp('Write in-plane file');
        obj.writeData(obj.inpFile,Minp);
        
        disp('Write out-of-plane file');
        obj.writeData(obj.outFile,Mout);
        
        obj.checkDeflection();
    end 
    
    % calculate components of demagnetizing tensor using in-plane and
    % out-of-plane components of magnetizationf and magnetic field
    function calcDemagTensor(obj,varargin)
        
        p = inputParser;
        p.addParamValue('field', 'effective', @(x) any(strcmp({'effective','demag','demag+exchange'},x)));
        p.parse(varargin{:});
        params = p.Results;
        
        % load static magnetization
        M0 = obj.convM*obj.readFile(obj.staticFile);
        Sz = size(M0);
        L = Sz(1)*Sz(2)*Sz(3);
        Ms = sqrt(M0(:,:,:,1).^2+M0(:,:,:,2).^2+M0(:,:,:,3).^2);
        obj.Ms = Ms;
        
        % load dynamic magnetization
        Minp = obj.convM*obj.readFile(obj.inpFile);
        Mout = obj.convM*obj.readFile(obj.outFile);
        
        switch (params.field)
            case 'effective'
                % load effective fields
                Heff0 = obj.convH*obj.readFile(obj.BdemagStaticFile);
                HeffInp = obj.convH*obj.readFile(obj.BexchStaticFile);
                HeffOut = obj.convH*obj.readFile(obj.BeffFile);
                
                % calculate Nzz
                obj.Nzz = -(Heff0(:,:,:,1).*M0(:,:,:,1)+Heff0(:,:,:,2).*M0(:,:,:,2)...
                    +Heff0(:,:,:,3).*M0(:,:,:,3))./(Ms.^2);
                
                % calculate dynamic fields
                hInp = reshape(HeffInp - Heff0,[L 3]);
                hOut = reshape(HeffOut - Heff0,[L 3]);
            
            case 'demag'
                % load demagnetizing fields
                Hd0   = obj.convH*obj.readFile(obj.BdemagStaticFile);
                HdInp = obj.convH*obj.readFile(obj.BdemagInpFile);
                HdOut = obj.convH*obj.readFile(obj.BdemagOutFile);
                
                % calculate Nzz
                obj.Nzz = -(Hd0(:,:,:,1).*M0(:,:,:,1)+Hd0(:,:,:,2).*M0(:,:,:,2)...
                    +Hd0(:,:,:,3).*M0(:,:,:,3))./(Ms.^2);
                
                % calculate dynamic fields
                hInp = reshape(HdInp - Hd0,[L 3]);
                hOut = reshape(HdOut - Hd0,[L 3]);
                
            case 'demag+exchange'
                % load demagnetizing fields
                Hd0   = obj.convH*obj.readFile(obj.BdemagStaticFile);
                HdInp = obj.convH*obj.readFile(obj.BdemagInpFile);
                HdOut = obj.convH*obj.readFile(obj.BdemagOutFile);
                
                % load exchange field
                Hexch0   = obj.convH*obj.readFile(obj.BexchStaticFile);
                HexchInp = obj.convH*obj.readFile(obj.BexchInpFile);
                HexchOut = obj.convH*obj.readFile(obj.BexchOutFile);
                
                % calculate Nzz
                obj.Nzz = -((Hd0(:,:,:,1)+Hexch0(:,:,:,1)).*M0(:,:,:,1)...
                    +(Hd0(:,:,:,2)+Hexch0(:,:,:,2)).*M0(:,:,:,2)...
                    +(Hd0(:,:,:,3)+Hexch0(:,:,:,3)).*M0(:,:,:,3))./(Ms.^2);
                
                % calculate dynamic fields
                hInp = reshape(HdInp - Hd0 + HexchInp - Hexch0,[L 3]);
                hOut = reshape(HdOut - Hd0 + HexchOut - Hexch0,[L 3]);
        end        

        
        mInp = reshape(Minp-M0,[L 3]);
        mOut = reshape(Mout-M0,[L 3]);
        M0   = reshape(M0,[L 3]);
        Ms   = reshape(Ms,[L 1]);
        
        parfor ind = 1:size(hInp,1)
            rot = obj.getRotation(M0(ind,1),M0(ind,2),M0(ind,3));
            
            % calculate dynamic magnetization in local system
            mInpLoc(ind,:) = rot*mInp(ind,:).';
            mOutLoc(ind,:) = rot*mOut(ind,:).';
            M0Loc(ind,:) = rot*M0(ind,:).';
                        
            % calculate dynamic fields in local coordinate system
            hInpLoc(ind,:)  = rot*hInp(ind,:).';
            hOutLoc(ind,:)  = rot*hOut(ind,:).';   
        end
        
        
        hOutLoc = reshape(hOutLoc,[Sz(1),Sz(2),Sz(3),3]);
        hInpLoc = reshape(hInpLoc,[Sz(1),Sz(2),Sz(3),3]);
        
        M0Loc = reshape(M0Loc,[Sz(1),Sz(2),Sz(3),3]);
        mOutLoc = reshape(mOutLoc,[Sz(1),Sz(2),Sz(3),3]);
        mInpLoc = reshape(mInpLoc,[Sz(1),Sz(2),Sz(3),3]);
                
        % calculate demagnetizing factors

        obj.Nyy = -hOutLoc(:,:,:,3)./mOutLoc(:,:,:,3);
        obj.Nyx = -hOutLoc(:,:,:,2)./mOutLoc(:,:,:,3);
        
        obj.Nxx = -hInpLoc(:,:,:,2)./mInpLoc(:,:,:,2);
        obj.Nxy = -hInpLoc(:,:,:,3)./mInpLoc(:,:,:,2);

        obj.calcFreq();
        obj.show();
                
        figure(10);
        imagesc(mInpLoc(:,:,1,2).');
        title('M inp Loc'); colorbar();
        
        figure(11);
        imagesc(hInpLoc(:,:,1,1).');
        title('H inp Loc X'); colorbar()
        
        figure(12);
        imagesc(hInpLoc(:,:,1,2).');
        title('H inp Loc Y'); colorbar();
        
        figure(13);
        imagesc(hInpLoc(:,:,1,3).');
        title('H inp Loc Z'); colorbar();

        
    end
    
    % calculate frequency of quasi-uniform mode using demagnetizing tensor 
    function calcFreq(obj)
        obj.freq = obj.gamma*sqrt(abs((obj.H+obj.Ms.*(obj.Nxx-obj.Nzz)).*...
            (obj.H+obj.Ms.*(obj.Nyy-obj.Nzz))))/(2*pi);
        
    end    
    
    % show results of calculation of demag tensor and frequency
    function show(obj,varargin)   
        % parse input parameters
        p = inputParser;
        p.addParamValue('saveImg', false,@islogical);
        p.addParamValue('zSlice',1,@isnumeric);
        p.addParamValue('saveMatAs','',@isstr);
        p.parse(varargin{:});
        params = p.Results;
        
        % check input parameters
        if params.zSlice > size(obj.Nxx,3)
            params.zSlice = 1;
            disp('Warning: parameter "zSlice" exceeds array size. It was set to 1.')
        end    
        
        % calculate scales
        xScale = linspace(obj.xmin,obj.xmax,obj.xnodes)*1e6;
        yScale = linspace(obj.ymin,obj.ymax,obj.ynodes)*1e6;

        
        f1=figure(1);
            clf(); set(gcf,'name','Nxx','numbertitle','off') 
            imagesc(xScale,yScale,squeeze(obj.Nxx(:,:,params.zSlice)).');
            title('Nxx'); colorbar();

        f2 = figure(2);
            clf(); set(gcf,'name','Nyy','numbertitle','off') 
            imagesc(xScale,yScale,squeeze(obj.Nyy(:,:,params.zSlice)).');
            title('Nyy'); colorbar();

        f3 = figure(3);
            clf(); set(gcf,'name','Nzz','numbertitle','off') 
            imagesc(xScale,yScale,squeeze(obj.Nzz(:,:,params.zSlice)).');
            title('Nzz'); colorbar();

        f4 = figure(4);
            clf(); set(gcf,'name','Frequency','numbertitle','off');
            imagesc(xScale,yScale,squeeze(obj.freq(:,:,params.zSlice)).'/1e9);
            title('Frequency'); colorbar();
        
        yScale = linspace(0,5,obj.ynodes);
        f5 = figure(5);
        clf();
        set(gcf,'name','Frequency slice','numbertitle','off'); 
        plot(yScale,squeeze(obj.freq(1000,:,1))/1e9); title('Frequency');
        ylabel('Frequency (GHz)');
        xlabel('x (\mum)');
        
        % save images
        if true
            print(f1,'-dpng','Nxx.png');
            print(f2,'-dpng','Nyy.png');
            print(f3,'-dpng','Nzz.png');
            print(f4,'-dpng','fMap.png');
            print(f5,'-dpng','fSlice.png');
        end
        
        % save frequency array to *.mat file
        if ~strcmp(params.saveMatAs,'')
            freq = obj.freq;
            save(strcat(params.saveMatAs,'.mat'),'freq');
        end    
    end   
    
    % plot result of deflection (test method)
    function checkDeflection(obj)
        slice = 1;
        % load static magnetization
        M0 = obj.convM*obj.readFile(obj.staticFile);
        Sz = size(M0);
        Ms = sqrt(M0(:,:,:,1).^2+M0(:,:,:,2).^2+M0(:,:,:,3).^2);
        figure(12)
        imagesc(Ms(:,:,1));
        
        % load dynamic magnetization
        Minp = obj.convM*obj.readFile(obj.inpFile);
        MsInp = sqrt(Minp(:,:,:,1).^2+Minp(:,:,:,2).^2+Minp(:,:,:,3).^2);
        figure(13)
        imagesc(MsInp(:,:,1));
        
        
        Mout = obj.convM*obj.readFile(obj.outFile);
        MsOut = sqrt(Mout(:,:,:,1).^2+Mout(:,:,:,2).^2+Mout(:,:,:,3).^2);
        figure(14)
        imagesc(MsOut(:,:,1));
             
        
        L = Sz(1)*Sz(2)*Sz(3);
        
        
        figure(6);
        clf();
        set(gcf,'name','M0','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(M0(:,slice,:,1)).');
            colorbar(); title('M_0^x');
        subplot(312)
            imagesc(squeeze(M0(:,slice,:,2)).');
            colorbar(); title('M_0^y')
        subplot(313)
            imagesc(squeeze(M0(:,slice,:,3)).');
            colorbar(); title('M_0^z');
            
        figure(7);
        clf();
        set(gcf,'name','Minp','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(Minp(:,slice,:,1)).');
            colorbar(); title('M_inp^x');
        subplot(312)
            imagesc(squeeze(Minp(:,slice,:,2)).');
            colorbar(); title('M_inp^y')
        subplot(313)
            imagesc(squeeze(Minp(:,slice,:,3)).');
            colorbar(); title('M_inp^z');
            
        figure(8);
        clf();
        set(gcf,'name','Mout','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(Mout(:,slice,:,1)).');
            colorbar(); title('M_out^x');
        subplot(312)
            imagesc(squeeze(Mout(:,slice,:,2)).');
            colorbar(); title('M_out^y')
        subplot(313)
            imagesc(squeeze(Mout(:,slice,:,3)).');
            colorbar(); title('M_out^z');    
        
        M0   = reshape(M0,[L 3]);
        Ms   = reshape(Ms,[L 1]);
        Minp = reshape(Minp,[L 3]);
        Mout = reshape(Mout,[L 3]);
        
        

        parfor ind = 1:size(Minp,1)
            rot = obj.getRotation(M0(ind,1),M0(ind,2),M0(ind,3));            
            % calculate magnetization in local system
            M0Loc(ind,:) = rot*M0(ind,:).';
            MinpLoc(ind,:) = rot*Minp(ind,:).';
            MoutLoc(ind,:) = rot*Mout(ind,:).';        
        end
        
        M0Loc = reshape(M0Loc,[Sz(1),Sz(2),Sz(3),3]);
        MinpLoc = reshape(MinpLoc,[Sz(1),Sz(2),Sz(3),3]);
        MoutLoc = reshape(MoutLoc,[Sz(1),Sz(2),Sz(3),3]);
        
        mInp = MinpLoc - M0Loc;
        mOut = MoutLoc - M0Loc;
        

        figure(1);
        clf();
        set(gcf,'name','M0 local','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(M0Loc(:,slice,:,1)).');
            colorbar(); title('M_0^z');
        subplot(312)
            imagesc(squeeze(M0Loc(:,slice,:,2)).');
            colorbar(); title('M_0^x')
        subplot(313)
            imagesc(squeeze(M0Loc(:,slice,:,3)).');
            colorbar(); title('M_0^y')
        
            
        figure(2);
        clf();
        set(gcf,'name','Minp Loc','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(MinpLoc(:,slice,:,1)).');
            colorbar(); title('M_{inp}^z');
        subplot(312)
            imagesc(squeeze(MinpLoc(:,slice,:,2)).');
            colorbar(); title('M_{inp}^x');
        subplot(313)
            imagesc(squeeze(MinpLoc(:,slice,:,3)).');
            colorbar(); title('M_{inp}^y');
            
        figure(3);
        clf();
        set(gcf,'name','Mout Loc','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(MoutLoc(:,slice,:,1)).');
            colorbar(); title('M_{out}^z');
        subplot(312)
            imagesc(squeeze(MoutLoc(:,slice,:,2)).');
            colorbar(); title('M_{out}^x');
        subplot(313)
            imagesc(squeeze(MoutLoc(:,slice,:,3)).');
            colorbar(); title('M_{out}^y');
            
        figure(4);
        clf();
        set(gcf,'name','mInp Loc','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(mInp(:,slice,:,1)).');
            colorbar(); title('M_{out}^x');
        subplot(312)
            imagesc(squeeze(mInp(:,slice,:,2)).');
            colorbar(); title('M_{out}^y');
        subplot(313)
            imagesc(squeeze(mInp(:,slice,:,3)).');
            colorbar(); title('M_{out}^z');
       
       figure(5);
        clf();
        set(gcf,'name','mOut Loc','numbertitle','off') 
        subplot(311)
            imagesc(squeeze(mOut(:,slice,:,1)).');
            colorbar(); title('M_{out}^x');
        subplot(312)
            imagesc(squeeze(mOut(:,slice,:,2)).');
            colorbar(); title('M_{out}^y');
        subplot(313)
            imagesc(squeeze(mOut(:,slice,:,3)).');
            colorbar(); title('M_{out}^z');     
        
    end
    
    % calculate spatial distribution of magnetic susceptibility for the
    % given frequency
    % parameters:
    %     - freq - frequency (GHz)
    function calcSuscept(obj,varargin)      
        % parse input parameters
        p = inputParser();
        p.addParamValue('freq',10,@isnumeric);
        p.addParamValue('showResults',true,@islogical);
        p.parse(varargin{:});
        params = p.Results;
        
        w = 2*pi*params.freq*1e9;
        
        wH = obj.gamma*(obj.H-obj.Nzz.*obj.Ms);
        w1 = wH +0.5*obj.gamma*obj.Ms.*(obj.Nxx+obj.Nyy+obj.Nyx+obj.Nxy);
        w0 = sqrt((obj.H+obj.Ms.*(obj.Nxx-obj.Nzz)).*(obj.H+obj.Ms.*(obj.Nyy-obj.Nzz))...
         -obj.Nxy.*obj.Nyx.*obj.Ms.^2).*obj.gamma;
     
        % (w0^2-w^2) denominator
        freqDiff = w0.^2-w.^2;
                
        chiXX = ((obj.gamma*obj.Ms.*freqDiff.*(wH+obj.gamma*obj.Ms.*obj.Nyy))...
            -j*(-obj.alpha*w*obj.gamma*obj.Ms.*freqDiff+...
            2*obj.alpha*obj.gamma*w*w1.*obj.Ms.*(wH+obj.gamma*obj.Ms.*obj.Nyy)))./freqDiff.^2;
        
        chiYY = ((obj.gamma*obj.Ms.*freqDiff.*(wH+obj.gamma*obj.Ms.*obj.Nxx))...
            -j*(-obj.alpha*w*obj.gamma*obj.Ms.*freqDiff+...
            2*obj.alpha*obj.gamma*w*w1.*obj.Ms.*(wH+obj.gamma*obj.Ms.*obj.Nxx)))./freqDiff.^2;
        
        chiXY = (-(obj.gamma*obj.Ms).^2.*freqDiff.*obj.Nxy+2*obj.gamma*obj.Ms.*w1.*w^2*obj.alpha-...
                    j*(-w*obj.gamma*obj.Ms.*freqDiff-8*pi*(obj.gamma*obj.Ms).^2.*w1.*obj.Nxy.*w*obj.alpha))./freqDiff.^2;     
       
        chiYX = (-(obj.gamma*obj.Ms).^2.*freqDiff.*obj.Nyx-2*obj.gamma*obj.Ms.*w1.*w^2*obj.alpha-...
                    j*(w*obj.gamma*obj.Ms.*freqDiff-8*pi*(obj.gamma*obj.Ms).^2.*w1.*obj.Nyx.*w*obj.alpha))./freqDiff.^2;     
        
        obj.suscept(params.freq) = struct('xx',chiXX,'xy',chiXY,'yx',chiYX,'yy',chiYY);
        
        if params.showResults
            % calculate scales
            xScale = linspace(obj.xmin,obj.xmax,obj.xnodes)*1e6;
            yScale = linspace(obj.ymin,obj.ymax,obj.ynodes)*1e6;
            
            figure(1);
                clf(); set(gcf,'name','chi_x_x','numbertitle','off')
                subplot(211);
                    imagesc(xScale,yScale,log10(abs(real(squeeze(obj.chiXX(:,:,1))))).');
                    colorbar()
                subplot(212);
                    imagesc(xScale,yScale,log10(abs(imag(squeeze(obj.chiXX(:,:,1))))).')
                    colorbar();

            figure(2);
                clf(); set(gcf,'name','chi_y_y','numbertitle','off')
                subplot(211);
                    imagesc(xScale,yScale,real(squeeze(obj.chiYY(:,:,1))).');
                    colorbar()
                subplot(212);
                    imagesc(xScale,yScale,imag(squeeze(obj.chiYY(:,:,1))).')
                    colorbar();

            figure(3);
                clf(); set(gcf,'name','chi_x_y','numbertitle','off')
                subplot(211);
                    imagesc(xScale,yScale,real(squeeze(obj.chiXY(:,:,1))).');
                    colorbar()
                subplot(212);
                    imagesc(xScale,yScale,imag(squeeze(obj.chiXY(:,:,1))).')
                    colorbar();    

            figure(4);
                clf(); set(gcf,'name','chi_y_x','numbertitle','off')
                subplot(211);
                    imagesc(xScale,yScale,real(squeeze(obj.chiYX(:,:,1))).');
                    colorbar()
                subplot(212);
                    imagesc(xScale,yScale,imag(squeeze(obj.chiYX(:,:,1))).')
                    colorbar();

            figure(5);
                hist(sqrt(abs(freqDiff(:))),20); title('w^2-w0^2');

            figure(6);
                hist(w0(:),20); title('w0');
        end         
    end
    
    % calculate set of spatial maps of magnetic susceptibility for the
    % given range of frequency
    function getSuscept(obj,freqList)
       for freq = freqList
           obj.calcSuscept('freq',freq);
       end    
    end    
    
end

methods (Access = private)
    
    % read data from *.ovf file
    % params:
    %     fName - name of input file
    function M = readFile(obj,fName)       
        fid = fopen(fName);
        [IOmess, errnum] = ferror(fid);
        if (errnum ~= 0)
            disp(IOmess);
            return;
        end
        
        % read parameters file
        expr = '^#\s([\w\s:]+):\s([-.0-9e]+)';
        propertiesList = fieldnames(obj);
        obj.header = '';
        line = fgetl(fid);
        obj.header = strcat(obj.header,line,'\n');
        while (isempty(strfind(line,'Begin: Data Binary')))
            line = fgetl(fid);
            obj.header = strcat(obj.header,line,'\n');
            [~, ~, ~, ~, tokenStr, ~, splitStr] = regexp(line,expr);
            % read parameters
            if (size(tokenStr,1)>0)
                if (size(tokenStr{1,1},2)>1)
                    % seek properties
                    toks = tokenStr{1,1};
                    
                    if (strcmp(toks{1,1},'Desc:  Iteration'))
                        %obj.iteration = str2num(toks{1,2});
                    elseif (strcmp(toks{1,1},'Desc:  Total simulation time'))
                        %obj.totalSimTime = str2num(toks{1,2});
                    else
                        for i=1:size(propertiesList,1)
                            if(strcmp(propertiesList{i,1},toks{1,1}))
                                prop = toks{1,1};
                                val = toks{1,2};
                                
                                %  Is it numerical value?
                                [num,status] = str2num(val);
                                if (status) % yes, it's numerical
                                    set(obj,prop,num)
                                else % no, it's string
                                    set(obj,prop,val)
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % determine file format

        if (~isempty(strfind(line,'8')))
            obj.format = 'double';
            testVal = obj.doubleTestVal;
        elseif (~isempty(strfind(line,'4')))
            obj.format = 'single';
            testVal = obj.singleTestVal;
        else
            disp('Unknown format');
            return
        end
        
        % read first test value
        fTestVal = fread(fid, 1, obj.format, 0, 'ieee-le');
        if ((strcmp(obj.format,'single') && (fTestVal == obj.singleTestVal))...
                || (strcmp(obj.format,'double') && (fTestVal == obj.doubleTestVal)))
            disp('Correct format')
        else
            disp('Wrong format');
            return;
        end
        
        
        data = fread(fid, obj.xnodes*obj.ynodes*obj.znodes*obj.dim,...
            obj.format, 0, 'ieee-le');
        
        line = fgetl(fid);
        if (isempty(strfind(line ,'# End: Data')) && isempty(strfind(line,'# End: Segment')))
            disp('End of file is incorrect. Something wrong');
            fclose(fid);
            % return;
        else
            fclose(fid);
        end

        Mx = data(1:3:size(data,1));
        My = data(2:3:size(data,1));
        Mz = data(3:3:size(data,1));
        
        M(:,:,:,1) = reshape(Mx, [obj.xnodes obj.ynodes obj.znodes]);
        M(:,:,:,2) = reshape(My, [obj.xnodes obj.ynodes obj.znodes]);
        M(:,:,:,3) = reshape(Mz, [obj.xnodes obj.ynodes obj.znodes]);
    end
    
    % write data to file
    % params:
    %    fName - name of output file
    %    arr   - array of data 
    function writeData(obj,fName,arr)
        data =[];
        for zInd = 1:obj.znodes
            for yInd = 1:obj.ynodes
                for xInd = 1:obj.xnodes
                    for projInd = 1:obj.dim
                        data = [data arr(xInd,yInd,zInd,projInd)];
                        %disp([num2str(xInd) ' ' num2str(yInd) ' ' num2str(zInd)])
                    end
                end
            end
        end
        
        try
            [fid,errMsg] = fopen(fName,'w');
            fprintf(fid,sprintf(obj.header));
            if strcmp(obj.format,'single')
                fwrite(fid, single(obj.singleTestVal), 'single', 0, 'ieee-le');
                fwrite(fid, single(data), 'single', 0, 'ieee-le');
                fprintf(fid,sprintf('\n#End: Data Binary 4\n#End: Segment\n'));
            elseif strcmp(obj.format,'double')
                fwrite(fid, obj.doubleTestVal, 'double', 0, 'ieee-le');
                fwrite(fid, data, 'double', 0, 'ieee-le');
                fprintf(fid,sprintf('\n#End: Data Binary 8\n#End: Segment\n'));                
            else
                disp('Unknown format')
            end    
            fclose(fid);
        catch err
            fclose(fid);
            disp(err.message);
        end
    end
    
    % return rotation matrix for one spatial point of the mesh
    function rot = getRotation(obj,Mx,My,Mz)
        Ms = sqrt(Mx^2+My^2+Mz^2);
        sq = sqrt(Mx^2+My^2);
        rot = [Mx/Ms, My/Ms, Mz/Ms;...
                 -My/sq, Mx/sq, 0;...
                 -Mx*Mz/(Ms*sq), -My*Mz/(Ms*sq), sq/Ms];
    end    
end    

end    