classdef Simulation
    %SIMULATION Phragmoplast Simulation Class
    %   This class will run the simulator and produce the figures and data

    properties
        Name
        SerialNumber
        Params
        Fractions
        Images
        ImageFiles
        Kymograph
        Turnovers
        ImageData
        Lengths
        RelativeTubulin
        Luminosity
        InitTable
        Init
        LuminosityTable
        TempFolder
        DataFolder
        HasRun;
        PlotNames = ["Plot_Fractions", "Plot_Lengths", "Plot_Seed_Locations", "Plot_Seed_Distribution", "Plot_Length_and_Tubulin", "Plot_Kymograph", "Plot_Turnovers"];

    end

    methods
        function this = Simulation(pars, name, sn)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            this.Name = name;
            this.HasRun = false;
            this.SerialNumber = sn;
            this.TempFolder = './temp/';
            this.DataFolder = append(this.TempFolder, '/DATA_', num2str(sn), '/');
            mkdir(this.TempFolder);
            mkdir(append(this.TempFolder, 'Params/'));
            mkdir(append(this.TempFolder, '/DATA_', num2str(sn), '/'));
            this.Params = {
                'fname_prefix', 'FrameNamePrefix','The data prefix',append(this.TempFolder, 'DATA_', num2str(sn), '/', name, '_');
                'DistType', 'DistributionType','The distribution type', 0; %2
                'NoMature_mode','No Mature Mode', 'Allows microtubules to enter a mature state',0; %3
                'Thickness', 'Thickness','Thickness of the phragmoplast in microns',1.5000; %4
                'Width', 'Width','Width of the phragmoplast in microns',5.0000; %5
                'dx', 'Dx','Size step in microns',0.02; %6
                'N_MT', 'NumMT','Number of microtubules',1000; %7
                'Npixel_X','Npixel_X','Number of X pixels in images',76*3; %8
                'Npixel_Y', 'Npixel_Y','Number of Y pixels in images',256*3; %9
                't_min','TMin','Minimum time of the simulation',-300; %10
                't_max','TMax','Maximum time of the simulation',300; %11
                'FRAP_t', 'FRAPTime','FRAP time',50; %12
                'FRAP_x1', 'FRAPX1','FRAP X1 coordinate',0.20; %13
                'FRAP_y1', 'FRAPY1','FRAP Y1 coordinate',0.50; %14
                'FRAP_x2', 'FRAPX2','FRAP X2 coordinate',1.30; %15
                'FRAP_y2', 'FRAPY2','FRAP Y2 coordinate',4.50; %16
                'bin_N', 'numBins','Number of bins in the histogram',25; %17
                'outputJPEG', 'OutputJPEG','Outputs jpeg files to the temporary folder (very slow)',0; %18
                'outputBMA', 'OutputBMA','Outputs bma image files to be read by matlab (fast)',1; %19
                'blurRadiusJPEG', 'BlurRadius','Image Blur Radius',1; %20
                'init_type', 'InitType','Type of initialization: RND or FIXED','RND'; %21
                'writeImageTime', 'NoFigureTime','Do not write the time on the figures',0; %22
                'rnd_seed', 'RandomSeed','Random seed to initiate simulation with',-1; %23
                'SimIndex', 'SimIndex', '', sn; %24

                'r_polym', 'RatePolymerization', 'Rate of polymerization in micron/s', 0.06; %25
                'r_depolym', 'RateDePolymerization', 'Rate of depolymerization in micron/s', 0.12; %26
                'r_sg', 'RateShrink2Grow','Rate of transition from shrinking to growing',0.059; %27
                'r_gs', 'RateGrow2Shrink','Rate of transition from growing to shrinking',0.012; %28
                'r_pg', 'RatePause2Grow','Rate of transition from pause to growing',0.084; %29
                'r_ps', 'RatePause2Shrink','Rate of transition from pause to shrinking',0.0115; %30
                'r_gp', 'RateGrow2Pause','Rate of transition from growth to pause',0.1141; %31
                'r_sp', 'RateShrink2Pause','Rate of transition from shrinking to pause',0.0154; %32
                'r_ps_CP', 'RatePause2ShrinkCP','Rate of transition from pause to shrinking in the defined cellplate region',0.013; %33
                'r_ps_DZ', 'RatePause2ShrinkDZ','Rate of transition from pause to shriknking in the defined distal region',1; %34

                'dens_MT', 'TubulinDensisty','Density of available tubulin per unit area',45.79; %35
                'L0', 'InitialLength','Maximum initial length or fixed length at seeding',0.805;  %36
                'r_reseed', 'RateReseed','Rate of reseeding',1; %37
                'frac_seed_middle_slope', 'MiddleSeedSlope','Slope of the seeds in the middle',-1; %38
                'frac_grow_out', 'FractionGrowOut','Fraction that grows out',0.05; %39
                'frac_grow_in', 'FractionGrowIn','Fraction that grows in',0.95; %40
                'frac_seed_distal', 'FractionSeedDistal','Fraction of seeds in the distal region',.0; %41
                'frac_seed_middle', 'FractionSeedMiddle','Fraction of seeds in the middle region',.0; %42
                'frac_seed_CP', 'FractionSeedCP','Fraction of seeds in the cell plate region',.0; %43
                'frac_treadmill', 'FractionTreadmill ','Fraction of seeds that treadmill (+ and - end dynamic)',0; %44
                'theta_max', 'ThetaMax','Maximum angle of microtubules',20; %45
                'CP_X_MAX', 'CellPlateXMAX','Maximum x value for seeds located in the cell plate region',0.05; %46
                'DZ_X_MIN', 'DistalXMIN','Mimimum x value for distal region effect',1.40; %47
                'DZ_SEED_MIN', 'DistalSeedMIN','Minimum x value for distal region seeds',1.5; %48
                'grf_dt', 'ImageTimeStep','Time step between images',2; %49
                'grf_max', 'ImageMaxLuminosity','Maximum image luminosity (scaled so G=1)',1; %50

                'kymograph_amp', 'KymographAmp','Kymograph luminosity amplitude',-0.8; %51

                'r_me_polym', 'RateMinusPolymerization','Rate of minus-end polymerization in micron/s',0.0327;
                'r_me_depolym', 'RateMinusDepolymerization','Rate of minus-end depolymerization in micron/s',0.046;
                'r_me_gs', 'RateMinusGrow2Shrink','Rate of minus-end transition from growing to shrinking',0.1;
                'r_me_sg', 'RateMinusShrink2Grow','Rate of minus-end transition from shrinking to growing',0.039;
                'r_me_pg', 'RateMinusPause2Grow','Rate of minus-end transition from pause to growing',0.0065;
                'r_me_ps', 'RateMinusPause2Shrink','Rate of minus-end transition from pause to shrinking',0.033;
                'r_me_gp', 'RateMinusGrow2Pause','Rate of minus-end transition from growing to pause',0.056;
                'r_me_sp', 'RateMinusShrink2Pause','Rate of minus-end transition from shrinking to pause',0.072;

                'TaxolCorr', 'TaxolCorrection','Taxol correction',1.0;
                'TdensCorr', 'DensityCorrection','Density correction parameter',1;

                'SeedParameterA', 'SeedMean','The mean of the seeding distribution',1.1;
                'SeedParameterB', 'SeedSpread','The standard deviation or spead of the seeding distributiob',.35;
                'reseed_slope','ReseedSlope','The slope of the reseed rate curve',-.66;
                'reseed_b', 'ReseedIntercept','The Y-Intercept of the reseed rate curve',1};

            if ~isempty(pars)
                this.Params(:,4) = pars;
            end

            this.Fractions = [];
            this.Images = {};
            this.ImageFiles = [];
            this.Turnovers = [];
            this.Lengths = [];
            this.RelativeTubulin = [];
            this.Luminosity = [];
            this.InitTable = [];
            this.Init = [];
            this.LuminosityTable = [];

        end


        function this = SetPars(pars)
            this.Params(:,4) = pars;
        end

        function [this] = GetPlot(this, name, parent)
            if this.HasRun == false
                return
            end
            if name == this.PlotNames(1)
                this = this.PlotFractions(parent);
            elseif name == this.PlotNames(2)
                this = this.PlotLengths(parent);
            elseif name == this.PlotNames(3)
                this = this.PlotSeedLocations(parent);
            elseif name == this.PlotNames(4)
                this = this.PlotSeedDistribution(parent);
            elseif name == this.PlotNames(5)
                this = this.PlotLengthAndTubulin(parent);
            elseif name == this.PlotNames(6)
                this = this.PlotKymograph(parent);
            elseif name == this.PlotNames(7)
                this = this.PlotTurnovers(parent);
            else
                ME = MException('Simulation:GetPlot:NameNotFound', 'Figure name not found');
                throw(ME);
            end
        end

        function this = Run(this)
            w = waitbar(0,'Running Simulation...');
            w.UserData = [0, 100];
            filename = append(this.TempFolder, 'Params/SimNumber_', num2str(this.SerialNumber), '.par');

            writecell(this.Params(:, [1 4]), append(filename, '.txt'), 'Delimiter', ' ')
            unix(append('mv ', append(filename, '.txt '), filename))
            fileName = append(this.TempFolder, "Params/SimNumber_", int2str(this.SerialNumber), ".par");

            waitbar(.05, w);

            mkdir(this.TempFolder, "DATA_" + num2str(this.SerialNumber));
            argument1 = append('./php_simulator_3s', ' ', fileName);
            unix(argument1);

            waitbar(.5, w);

            fracFile = dir(fullfile(this.DataFolder, '*_Fractions.csv'));
            this.Fractions = readmatrix(fullfile(fracFile.folder, fracFile.name));

            waitbar(.55, w);

            kymoFile = dir(fullfile(this.DataFolder, '*_kymograph.jpeg'));
            this.Kymograph = imread(fullfile(kymoFile.folder, kymoFile.name));

            waitbar(.6, w);

            lenFile = dir(fullfile(this.DataFolder, '*_LengthsFile.csv'));
            this.Lengths = readmatrix(fullfile(lenFile.folder, lenFile.name), 'Delimiter', ',');

            waitbar(.65, w);

            this.ImageFiles = dir(fullfile(this.DataFolder, '*.bma'));
            this = this.ReadAllImages();

            waitbar(.80, w);

            tubFile = dir(fullfile(this.DataFolder, '*_RelativeMTDensity.csv'));
            this.RelativeTubulin = readmatrix(fullfile(tubFile.folder, tubFile.name), 'Delimiter', ',');

            waitbar(.85, w);

            initFile = dir(fullfile(this.DataFolder, '*InitFile.csv'));
            this.InitTable = readtable(fullfile(initFile.folder, initFile.name), 'Delimiter', ',', 'ReadVariableNames', true);
            this.Init = readmatrix(fullfile(initFile.folder, initFile.name), 'Delimiter', ',');


            gOutIdx = this.InitTable.type == 3;
            gInIdx = this.InitTable.type == 4;
            newData = {};
            newData(gOutIdx) = {'Grow-Out'};
            newData(gInIdx) = {'Grow-In'};
            this.InitTable.type = newData';

            waitbar(.90, w);

            lumFile = dir(fullfile(this.DataFolder, '*luminosity.txt'));
            this.LuminosityTable = readtable(fullfile(lumFile.folder, lumFile.name), 'Delimiter', ',', 'ReadVariableNames', true);
            this.Luminosity = readmatrix(fullfile(lumFile.folder, lumFile.name), 'Delimiter', ',');
            this.HasRun = true;

            waitbar(1, w);
            delete(w);
        end

        function this = PlotFractions(this, parent)


            averageFractions = mean(this.Fractions(:,4:end-1),1);


            %             histfit(EMTL, 25, 'Gamma');

            try
                parent.NextPlot = 'replace';
            catch
                hold off
            end

            histogram(parent, 'BinEdges', 0:.06:1.5, 'BinCounts', averageFractions);
            title(parent, 'Average of All Fractions');

        end

        function this = PlotLengths(this, parent)
            bins = [];
            lengths = this.Lengths;
            if ~isempty(gcp("nocreate"))
                parfor i = 1:size(this.Lengths, 1)
                    f = figure('visible','off');
                    hg = histogram(nonzeros(lengths(i,2:end)), 25, 'visible', 'off');
                    bins(i,:) = hg.Values;
                    delete(f);
                end
            else
                for i = 1:size(this.Lengths, 1)
                    f = figure('visible','off');
                    hg = histogram(nonzeros(lengths(i,2:end)), 25, 'visible', 'off');
                    bins(i,:) = hg.Values;
                    delete(f);
                end
            end

            try
                parent.NextPlot = 'replace';
            catch
                hold off
            end
            bar3(parent, bins)

        end

        function this = PlotKymograph(this, parent)
            parent.NextPlot = 'replace';
            imshow(this.Kymograph, 'Parent', parent)
        end

        function this = PlotSeedLocations(this, parent)
            grow_out = this.Init(this.Init(:,6) == 3, :);
            grow_in = this.Init(this.Init(:,6) == 4, :);

            try
                parent.NextPlot = 'replace';
            catch
                hold off
            end
            scatter(parent, grow_out(:,2), grow_out(:,3))
            try
                parent.NextPlot = 'add';
            catch
                hold on
            end
            scatter(parent, grow_in(:,2), grow_in(:,3))
            title(parent, 'MT Init Locations')
            legend(parent, 'Grow Out', 'Grow In', "Location", "bestoutside")

        end

        function this = PlotSeedDistribution(this, parent)
            try
                parent.NextPlot = 'replace';
            catch
                hold off
            end
            histogram(parent, this.Init(:,2))
            title(parent, 'X Seed Locations')
        end

        function this = PlotLengthAndTubulin(this, parent)

            try
                parent.NextPlot = 'replace';
            catch
                hold off
            end
            plot(parent, this.Luminosity(:,1), normalize(this.Luminosity(:,11)))
            try
                parent.NextPlot = 'add';
            catch
                hold on
            end
            plot(parent, this.Luminosity(:,1), normalize(this.Luminosity(:,13)))
            legend(parent, 'Total Length', 'RelativeMT')
        end



        function this = ReadAllImages(this)
            this.Images = struct();

            x1 = this.Params{13,4} * this.Params{8,4}/this.Params{4,4};
            y1 = this.Params{14,4} * this.Params{8,4}/this.Params{4,4};
            x2 = this.Params{15,4} * this.Params{8,4}/this.Params{4,4};
            y2 = this.Params{16,4} * this.Params{8,4}/this.Params{4,4};;

            blurRadius = this.Params{20, 4};
            imageFiles = this.ImageFiles;
            timeOnImage = this.Params{22,4};
            yA = [];
            yP = [];
            yD = [];
            if ~isempty(gcp("nocreate"))

                parfor i = 1:numel(this.ImageFiles)

                    file = fopen(fullfile(imageFiles(i).folder, imageFiles(i).name), 'r');
                    numX = fread(file, 1, 'int32');
                    numY = fread(file, 1, 'int32');
                    lenX = fread(file, 1, 'float64');
                    lenY = fread(file, 1, 'float64');
                    time = fread(file, 1, 'float64');
                    data = fread(file, numX * numY, 'float64');
                    fclose(file);
                    vmax = max(data);
                    data = data ./ vmax;

                    RGBArray = zeros(numY, numX, 3);
                    RGBArray(:,:,2) = reshape(data, numX, numY)';

                    % Blur the image if necessary
                    if blurRadius > 0
                        RGBArray = imgaussfilt(RGBArray, blurRadius);
                    end



                    Image = rgb2gray(RGBArray);

                    ImageCrop = imcrop(Image,[x1, y1, x2-x1, y2-y1]);
                    yA(i) = mean2(ImageCrop);

                    ImageCrop = imcrop(Image,[x1, y1, (x2-x1)/2, y2-y1]);
                    yP(i) = mean2(ImageCrop);

                    ImageCrop = imcrop(Image,[((x2 - x1)/2 + x1), y1, (x2-x1)/2, y2-y1]);
                    yD(i) = mean2(ImageCrop);




                    % Insert the image text
                    if timeOnImage == 0
                        RGBArray = insertText(RGBArray,[0, 0], string(time),'FontSize',18, 'Font', 'Arial', 'TextColor','red', 'BoxOpacity', 0, 'AnchorPoint', 'LeftTop');
                    end

                    images(i).RGBArray = RGBArray;
                    images(i).Time = time;

                end
                this.Images = images;
            else
                for i = 1:numel(this.ImageFiles)

                    file = fopen(fullfile(this.ImageFiles(i).folder, this.ImageFiles(i).name), 'r');
                    numX = fread(file, 1, 'int32');
                    numY = fread(file, 1, 'int32');
                    lenX = fread(file, 1, 'float64');
                    lenY = fread(file, 1, 'float64');
                    time = fread(file, 1, 'float64');
                    data = fread(file, numX * numY, 'float64');
                    fclose(file);
                    vmax = max(data);
                    data = data ./ vmax;

                    RGBArray = zeros(numY, numX, 3);
                    RGBArray(:,:,2) = reshape(data, numX, numY)';

                    % Blur the image if necessary
                    if blurRadius > 0
                        RGBArray = imgaussfilt(RGBArray, blurRadius);
                    end



                    Image = rgb2gray(RGBArray);

                    ImageCrop = imcrop(Image,[x1, y1, x2-x1, y2-y1]);
                    yA(i) = mean2(ImageCrop);

                    ImageCrop = imcrop(Image,[x1, y1, (x2-x1)/2, y2-y1]);
                    yP(i) = mean2(ImageCrop);

                    ImageCrop = imcrop(Image,[((x2 - x1)/2 + x1), y1, (x2-x1)/2, y2-y1]);
                    yD(i) = mean2(ImageCrop);




                    % Insert the image text
                    RGBArray = insertText(RGBArray,[0, 0], string(time),'FontSize',18, 'Font', 'Arial', 'TextColor','red', 'BoxOpacity', 0, 'AnchorPoint', 'LeftTop');

                    this.Images(i).RGBArray = RGBArray;
                    this.Images(i).Time = time;

                end
            end
            disp('Images loaded');
            [~,idx] = sort([this.Images.Time]);
            this.Images = this.Images(idx);

            yA = yA(idx);
            yP = yP(idx);
            yD = yD(idx);

            i = numel(this.ImageFiles);

            D = transpose(yD(26:i));
            P = transpose(yP(26:i));
            WR = transpose(yA(26:i));

            t = 0:2:i*2 + 2;
            t = transpose(t(26:i));
            t = t - t(1);

            ft = fittype( 'f_inf - (f_inf-f0)*exp(-k*x)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.MaxFunEvals = 900;
            opts.MaxIter = 500;
            opts.DiffMaxChange = 0.0001;
            opts.Lower = [min(D)-5 max(D)-5 .01];
            opts.Upper = [min(D)+5 max(D)+5 .1];

            [fitresult_D, ~] = fit( t, D, ft, opts );
            %
            opts.Lower = [min(P)-5 max(P)-5 .01];
            opts.Upper = [min(P)+5 max(P)+5 .1];

            [fitresult_P, ~] = fit( t, P, ft, opts );
            %
            opts.Lower = [min(WR)-5 max(WR)-5 .01];
            opts.Upper = [min(WR)+5 max(WR)+5 .1];

            [fitresult_WR, ~] = fit( t, WR, ft, opts);

            fit_values{1} = coeffvalues(fitresult_D);
            fit_values{2} = coeffvalues(fitresult_P);
            fit_values{3} = coeffvalues(fitresult_WR);
            K_Off = [fit_values{1}(3), fit_values{2}(3), fit_values{3}(3)];

            this.ImageData.Time = t;
            this.ImageData.D = D;
            this.ImageData.P = P;
            this.ImageData.WR = WR;
            this.ImageData.FitD = fitresult_D;
            this.ImageData.FitP = fitresult_P;
            this.ImageData.FitWR = fitresult_WR;
            this.Turnovers = log(2)./K_Off;


            % figure();

        end

        function this = PlotTurnovers(this, parent)

            t = this.ImageData.Time;
            D = this.ImageData.D;
            P = this.ImageData.P;
            WR = this.ImageData.WR;

            fitresult_D = this.ImageData.FitD;
            fitresult_P = this.ImageData.FitP;
            fitresult_WR = this.ImageData.FitWR;

            plot(parent, t, D);
            hold(parent);
            plot(parent, t, P)
            plot(parent, t, WR)
            plot(parent, t, fitresult_D(t));
            plot(parent, t,fitresult_P(t));
            plot(parent, t,fitresult_WR(t));
            legend(parent, 'Distal', 'Proximal',  'Whole', 'Fit_D', 'Fit_P','Fit_W');
            xlabel(parent, 'Time (s)')
            ylabel(parent, 'Luminosity')
        end

        function ExportMatFile(this, directory)
            save(append(directory, '/', this.Name, '_Package.mat'), "this",'-mat');
        end

        function ExportDataPackage(this, fileDir)
            if ~isfolder(fileDir)
                mkdir(fileDir)
            end
            fileDir = append(fileDir,'/',this.Name);
            mkdir(fileDir);

            mkdir(append(fileDir, '/Images'))
            for i = 1:numel(this.Images)
                imwrite(this.Images(i).RGBArray, append(fileDir, '/Images/', 'Image_', num2str(i), '.jpg'))
            end

            imwrite(this.Kymograph, append(fileDir, '/', 'Kymograph.jpeg'))

            writetable(this.InitTable, append(fileDir, '/InitTable.csv'))
            writetable(this.LuminosityTable, append(fileDir, '/LuminosityTable.csv'))
            writecell(this.Params, append(fileDir, '/Params.csv'))
            Info{1} = this.Name;
            Info{2} = this.SerialNumber;
            Info{3} = this.Turnovers;
            writecell(Info, append(fileDir, '/Info.csv'), "FileType", "text", "Delimiter",',');
            this.ExportAllPlots(fileDir)
        end

        function ExportPlot(this, fileDir, plotIndex)
            u = figure('visible', 'off');
            f = axes(u);
            this.GetPlot(plotIndex, f);
            saveas(u, append(fileDir, '/', this.Name, '_', plotIndex, '.jpeg'));
            delete(f)
            close(u)
        end

        function ExportAllPlots(this, fileDir)
            u = figure('visible', 'off');
            for i = 1:numel(this.PlotNames)
                f = axes(u);
                this.GetPlot(this.PlotNames{i}, f);
                saveas(u, append(fileDir, '/', this.Name, '_', this.PlotNames{i}, '.jpeg'));
                delete(f)
            end
            close(u)
        end

        function ExportParameters(this, fileDir)
            writecell(this.Params, append(fileDir, '/', this.Name, '_Params.csv'));
        end
    end
end

