classdef SimulatorApp_Export < matlab.apps.AppBase

    properties (Constant)
        HelpGuidePath = './Help/SimulatorHelp.pdf';
        DetailsInformation = 'This app simulates microtubule dynamics in the phragmoplast. It was developed by Bernard Piette (University of Durham?) and Matthew Hickey (Washington State University).Information on use, the code, and the parameters can be found in the help directory (./Help/)';
    end
    properties (Access = private)

        ParPool;            % Parallel Pool for parfor loops
        MaxThreads;         % Max thread count in pool
        Params;             % Current Params
        SimNames;           % List of simulation names
        Sim;                % List of simulations
        CurrentSim;         % Current sim index
        CurrentImageIndex;  % Current image index
        CurrentPlotName;    % Current plot index
        TableData;          % Main table data
        NumSims;            % Total number of sims
        CurrentDataTable;   % The current data table to display
        ImageTableData;     % The image table data

    end
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        DetailsMenu                     matlab.ui.container.Menu
        ExitMenu                        matlab.ui.container.Menu
        LoadMenu                        matlab.ui.container.Menu
        LoadSimulatorSessionmatMenu     matlab.ui.container.Menu
        LoadParametersFileMenu          matlab.ui.container.Menu
        LoadDataFromFolderMenu          matlab.ui.container.Menu

        ExportMenu                      matlab.ui.container.Menu
        ExportSimulatorSessionmatMenu   matlab.ui.container.Menu
        ExportParameterFilesMenu        matlab.ui.container.Menu
        ExportCurrentPlotMenu           matlab.ui.container.Menu
        ExportAllPlotsMenu              matlab.ui.container.Menu
        ExportRawDatatoFolderMenu       matlab.ui.container.Menu
        SimulationMenu                  matlab.ui.container.Menu
        RunMenu                         matlab.ui.container.Menu
        RunSettingsMenu                 matlab.ui.container.Menu
        MaxMemoryMenu                   matlab.ui.container.Menu
        MaxThreadsMenu                  matlab.ui.container.Menu
        TemporaryFileLocationMenu       matlab.ui.container.Menu
        CurrentMenu                     matlab.ui.container.Menu
        ChangeRandomSeedMenu            matlab.ui.container.Menu
        TabGroup                        matlab.ui.container.TabGroup
        ParametersTab                   matlab.ui.container.Tab
        InitTypeDropDown                matlab.ui.control.DropDown
        InitTypeDropDownLabel           matlab.ui.control.Label
        ImageBlurRadiusEditFieldLabel   matlab.ui.control.Label
        ImageBlurRadiusEditField        matlab.ui.control.NumericEditField
        SessionSimulationsDropDown      matlab.ui.control.DropDown
        SessionSimulationsDropDownLabel  matlab.ui.control.Label
        NumberofBinsEditFieldLabel      matlab.ui.control.Label
        NumberofBinsEditField           matlab.ui.control.NumericEditField
        FRAPY2EditFieldLabel            matlab.ui.control.Label
        FRAPY2EditField                 matlab.ui.control.NumericEditField
        FRAPX2EditFieldLabel            matlab.ui.control.Label
        FRAPX2EditField                 matlab.ui.control.NumericEditField
        FRAPY1EditFieldLabel            matlab.ui.control.Label
        FRAPY1EditField                 matlab.ui.control.NumericEditField
        FRAPX1EditFieldLabel            matlab.ui.control.Label
        FRAPX1EditField                 matlab.ui.control.NumericEditField
        FRAPTimeEditFieldLabel          matlab.ui.control.Label
        FRAPTimeEditField               matlab.ui.control.NumericEditField
        MaximumTimeEditFieldLabel       matlab.ui.control.Label
        MaximumTimeEditField            matlab.ui.control.NumericEditField
        MinimumTimeEditFieldLabel       matlab.ui.control.Label
        MinimumTimeEditField            matlab.ui.control.NumericEditField
        NumberofYPixelsEditFieldLabel   matlab.ui.control.Label
        NumberofYPixelsEditField        matlab.ui.control.NumericEditField
        NumberofXPixelsEditFieldLabel   matlab.ui.control.Label
        NumberofXPixelsEditField        matlab.ui.control.NumericEditField
        NumberofMTEditFieldLabel        matlab.ui.control.Label
        NumberofMTEditField             matlab.ui.control.NumericEditField
        ChangeinXEditFieldLabel         matlab.ui.control.Label
        ChangeinXEditField              matlab.ui.control.NumericEditField
        WidthEditFieldLabel             matlab.ui.control.Label
        WidthEditField                  matlab.ui.control.NumericEditField
        ThicknessEditFieldLabel         matlab.ui.control.Label
        ThicknessEditField              matlab.ui.control.NumericEditField
        DistributionTypeDropDownLabel   matlab.ui.control.Label
        DistributionTypeDropDown        matlab.ui.control.DropDown
        SimulationNameEditFieldLabel    matlab.ui.control.Label
        SimulationNameEditField         matlab.ui.control.EditField
        NoTimeonImagesButton            matlab.ui.control.StateButton
        OutputJPEGsButton               matlab.ui.control.StateButton
        OutputBMAsButton                matlab.ui.control.StateButton
        RunButton                       matlab.ui.control.Button
        ParameterTable                  matlab.ui.control.Table
        NoMatureModeButton              matlab.ui.control.StateButton
        ImagesTab                       matlab.ui.container.Tab
        NextPicture                     matlab.ui.control.Button
        PreviousPicture                 matlab.ui.control.Button
        NumPhotos                       matlab.ui.control.Label
        SessionSimulationsDropDown_2    matlab.ui.control.DropDown
        SessionSimulationsDropDown_2Label  matlab.ui.control.Label
        SimulationNameEditField_2       matlab.ui.control.EditField
        SimulationNameEditField_2Label  matlab.ui.control.Label
        UITable                         matlab.ui.control.Table
        Image                           matlab.ui.control.Image
        PlotsTab                        matlab.ui.container.Tab
        SessionSimulationsDropDown_3    matlab.ui.control.DropDown
        SessionSimulationsDropDown_3Label  matlab.ui.control.Label
        SimulationNameEditField_3       matlab.ui.control.EditField
        SimulationNameEditField_3Label  matlab.ui.control.Label
        PlotDropDown                    matlab.ui.control.DropDown
        PlotDropDownLabel               matlab.ui.control.Label
        DataDropDown                    matlab.ui.control.DropDown
        DataDropDownLabel               matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
        DataTab                         matlab.ui.container.Tab
        SessionSimulationsDropDown_4    matlab.ui.control.DropDown
        SessionSimulationsDropDown_4Label  matlab.ui.control.Label
        SimulationNameEditField_4       matlab.ui.control.EditField
        SimulationNameEditField_4Label  matlab.ui.control.Label
        UITable2                        matlab.ui.control.Table
    end

    % Callbacks that handle component events
    methods (Access = private)

        function UpdateImageGUI(app)
            if app.Sim{app.CurrentSim}.HasRun == true
                im = app.Sim{app.CurrentSim}.Images(app.CurrentImageIndex);
                app.Image.ImageSource = im.RGBArray;
                app.NumPhotos.Text = sprintf('%d of %d', app.CurrentImageIndex, numel(app.Sim{app.CurrentSim}.Images));

                app.UITable.Data = app.ImageTableData;
            end
        end

        function UpdateParameterGUI(app)
            app.NoMatureModeButton.Value = app.Sim{app.CurrentSim}.Params{3,4};
            app.ParameterTable.Data = app.Sim{app.CurrentSim}.Params(25:end, 2:4);
            app.WidthEditField.Value = app.Sim{app.CurrentSim}.Params{5,4};
            app.ThicknessEditField.Value = app.Sim{app.CurrentSim}.Params{4,4};
            app.ChangeinXEditField.Value = app.Sim{app.CurrentSim}.Params{6,4};
            app.NumberofMTEditField.Value = app.Sim{app.CurrentSim}.Params{7,4};
            app.NumberofXPixelsEditField.Value = app.Sim{app.CurrentSim}.Params{8,4};
            app.NumberofYPixelsEditField.Value = app.Sim{app.CurrentSim}.Params{9,4};
            app.MinimumTimeEditField.Value = app.Sim{app.CurrentSim}.Params{10,4};
            app.MaximumTimeEditField.Value = app.Sim{app.CurrentSim}.Params{11,4};
            app.FRAPTimeEditField.Value = app.Sim{app.CurrentSim}.Params{12,4};
            app.FRAPX1EditField.Value = app.Sim{app.CurrentSim}.Params{13,4};
            app.FRAPY1EditField.Value = app.Sim{app.CurrentSim}.Params{14,4};
            app.FRAPX2EditField.Value = app.Sim{app.CurrentSim}.Params{15,4};
            app.FRAPY2EditField.Value = app.Sim{app.CurrentSim}.Params{16,4};
            app.NumberofBinsEditField.Value = app.Sim{app.CurrentSim}.Params{17,4};
            app.OutputJPEGsButton.Value = app.Sim{app.CurrentSim}.Params{18,4};
            app.OutputBMAsButton.Value = app.Sim{app.CurrentSim}.Params{19,4};
            app.ImageBlurRadiusEditField.Value = app.Sim{app.CurrentSim}.Params{20,4};
            app.InitTypeDropDown.Value = app.Sim{app.CurrentSim}.Params{21,4};
            app.NoTimeonImagesButton.Value = app.Sim{app.CurrentSim}.Params{22,4};
        end

        function UpdateSimNameDropDowns(app)
            app.SessionSimulationsDropDown.Items = app.SimNames;
            app.SessionSimulationsDropDown_2.Items = app.SimNames;
            app.SessionSimulationsDropDown_3.Items = app.SimNames;
            app.SessionSimulationsDropDown_4.Items = app.SimNames;

            app.SessionSimulationsDropDown.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown_2.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown_3.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown_4.Value = app.SimNames(app.CurrentSim);

        end

        function UpdatePlotsGUI(app)
            if app.Sim{app.CurrentSim}.HasRun == true
                app.PlotDropDown.Items = app.Sim{app.CurrentSim}.PlotNames;
                app.Sim{app.CurrentSim}.GetPlot(app.CurrentPlotName, app.UIAxes);
            end
        end

        function DataDropDownValueChanged(app, ~)
            value = app.DataDropDown.Value;
            if strcmp(value, 'Turnovers') == 1
                app.CurrentDataTable = 1;
            elseif strcmp(value, 'Luminosity') == 1
                app.CurrentDataTable = 2;
            else
                app.CurrentDataTable = 3;
            end
            UpdateDataTable(app);
        end


        function UpdateDataTable(app)
            if app.CurrentDataTable == 1
                app.UITable2.Data = app.TableData;
                app.UITable2.ColumnName = {'Name', 'Whole (s)', 'Proximal (s)', 'Distal (s)'};
            elseif app.CurrentDataTable == 2
                app.UITable2.Data = app.Sim{app.CurrentSim}.LuminosityTable;
                app.UITable2.ColumnName = {'Time', 'Luminosity', '# Empty', '# Mature', '# Growing', '# Shrinking', '# Pause', '# Grow-Out', '# Grow-In', '# Treadmill', 'Total Length', 'Average Length', 'Relative Tubulin'};
            else
                app.UITable2.Data = app.Sim{app.CurrentSim}.InitTable;
                app.UITable2.ColumnName =  {'MT Index', 'Seed X-Location', 'Seed Y-Location', 'Initial Length', 'Theta', 'Type', 'Location'};
            end

        end

        function UpdateSimNameFields(app, value)
            app.Sim{app.CurrentSim}.Name = value;
            app.SimNames{app.CurrentSim} = value;
            app.UpdateSimNameDropDowns();

            app.SimulationNameEditField.Value = value;
            app.SimulationNameEditField_2.Value = value;
            app.SimulationNameEditField_3.Value = value;
            app.SimulationNameEditField_4.Value = value;

            app.Sim{app.CurrentSim}.Name = value;
            app.Sim{app.CurrentSim}.Params(1,4) = {append(app.Sim{app.CurrentSim}.TempFolder, 'DATA_', num2str(app.Sim{app.CurrentSim}.SerialNumber), '/', value, '_')};
        end

        % Value changed function: SimulationNameEditField
        function SimulationNameEditFieldValueChanged(app, ~)
            value = app.SimulationNameEditField.Value;
            app.UpdateSimNameFields(value);
        end

        % Value changed function: DistributionTypeDropDown
        function DistributionTypeDropDownValueChanged(app, ~)
            value = app.DistributionTypeDropDown.Value;
            if strcmp(value, 'Cauchy') == 1
                app.Sim{app.CurrentSim}.Params{2,4} = 0;
            elseif strcmp(value, 'Normal') == 1
                app.Sim{app.CurrentSim}.Params{2,4} = 1;
            elseif strcmp(value, 'Weibull') == 1
                app.Sim{app.CurrentSim}.Params{2,4} = 2;
            elseif strcmp(value, 'Exponential') == 1
                app.Sim{app.CurrentSim}.Params{2,4} = 3;
            else
                ME = MException('Simulator::DistributionDropdownChanged', 'Value out of range');
                throw(ME);
            end
        end

        % Value changed function: NoMatureModeButton
        function NoMatureModeButtonValueChanged(app, ~)
            value = app.NoMatureModeButton.Value;
            app.Sim{app.CurrentSim}.Params{3,4} = value;
        end

        % Cell edit callback: ParameterTable
        function ParameterTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            app.Sim{app.CurrentSim}.Params(indices(1)+24, 4) = {newData};
            app.ParameterTable.Data = app.Sim{app.CurrentSim}.Params(25:end, 2:4);
        end

        % Value changed function: WidthEditField
        function WidthEditFieldValueChanged(app, ~)
            value = app.WidthEditField.Value;
            app.Sim{app.CurrentSim}.Params{5,4} = value;
        end

        % Value changed function: ThicknessEditField
        function ThicknessEditFieldValueChanged(app, ~)
            value = app.ThicknessEditField.Value;
            app.Sim{app.CurrentSim}.Params{4,4} = value;
        end

        % Value changed function: ChangeinXEditField
        function ChangeinXEditFieldValueChanged(app, ~)
            value = app.ChangeinXEditField.Value;
            app.Sim{app.CurrentSim}.Params{6,4} = value;
        end

        % Value changed function: NumberofMTEditField
        function NumberofMTEditFieldValueChanged(app, ~)
            value = app.NumberofMTEditField.Value;
            app.Sim{app.CurrentSim}.Params{7,4} = value;
        end

        % Value changed function: NumberofXPixelsEditField
        function NumberofXPixelsEditFieldValueChanged(app, ~)
            value = app.NumberofXPixelsEditField.Value;
            app.Sim{app.CurrentSim}.Params{8,4} = value;
        end

        % Value changed function: NumberofYPixelsEditField
        function NumberofYPixelsEditFieldValueChanged(app, ~)
            value = app.NumberofYPixelsEditField.Value;
            app.Sim{app.CurrentSim}.Params{9,4} = value;
        end

        % Value changed function: MinimumTimeEditField
        function MinimumTimeEditFieldValueChanged(app, ~)
            value = app.MinimumTimeEditField.Value;
            app.Sim{app.CurrentSim}.Params{10,4} = value;
        end

        % Value changed function: MaximumTimeEditField
        function MaximumTimeEditFieldValueChanged(app, ~)
            value = app.MaximumTimeEditField.Value;
            app.Sim{app.CurrentSim}.Params{11,4} = value;
        end

        % Value changed function: FRAPTimeEditField
        function FRAPTimeEditFieldValueChanged(app, ~)
            value = app.FRAPTimeEditField.Value;
            app.Sim{app.CurrentSim}.Params{12,4} = value;
        end

        % Value changed function: FRAPX1EditField
        function FRAPX1EditFieldValueChanged(app, ~)
            value = app.FRAPX1EditField.Value;
            app.Sim{app.CurrentSim}.Params{13,4} = value;
        end

        % Value changed function: FRAPY1EditField
        function FRAPY1EditFieldValueChanged(app, ~)
            value = app.FRAPY1EditField.Value;
            app.Sim{app.CurrentSim}.Params{14,4} = value;
        end

        % Value changed function: FRAPX2EditField
        function FRAPX2EditFieldValueChanged(app, ~)
            value = app.FRAPX2EditField.Value;
            app.Sim{app.CurrentSim}.Params{15,4} = value;
        end

        % Value changed function: FRAPY2EditField
        function FRAPY2EditFieldValueChanged(app, ~)
            value = app.FRAPY2EditField.Value;
            app.Sim{app.CurrentSim}.Params{16,4} = value;
        end

        % Value changed function: NumberofBinsEditField
        function NumberofBinsEditFieldValueChanged(app, ~)
            value = app.NumberofBinsEditField.Value;
            app.Sim{app.CurrentSim}.Params{17,4} = value;
        end

        % Value changed function: OutputJPEGsButton
        function OutputJPEGsButtonValueChanged(app, ~)
            value = app.OutputJPEGsButton.Value;
            app.Sim{app.CurrentSim}.Params{18,4} = value;
        end

        % Value changed function: OutputBMAsButton
        function OutputBMAsButtonValueChanged(app, ~)
            value = app.OutputBMAsButton.Value;
            app.Sim{app.CurrentSim}.Params{19,4} = value;
        end

        % Value changed function: ImageBlurRadiusEditField
        function ImageBlurRadiusEditFieldValueChanged(app, ~)
            value = app.ImageBlurRadiusEditField.Value;
            app.Sim{app.CurrentSim}.Params{20,4} = value;
        end

        % Value changed function: InitTypeDropDown
        function InitTypeDropDownValueChanged(app, ~)
            value = app.InitTypeDropDown.Value;
            app.Sim{app.CurrentSim}.Params{21,4} = value;
        end

        % Value changed function: NoTimeonImagesButton
        function NoTimeonImagesButtonValueChanged(app, ~)
            value = app.NoTimeonImagesButton.Value;
            app.Sim{app.CurrentSim}.Params{22,4} = value;
        end


        % Menu selected function: ChangeRandomSeedMenu
        function ChangeRandomSeedMenuSelected(app, ~)
            input = inputdlg('Input the new random seed. Must be a positive integer or use -1 for the current time');
            app.Sim{app.CurrentSim}.Params{23,4} = input;
        end

        % Menu selected function: HelpMenu
        function HelpMenuSelected(app, ~)
            open(app.HelpGuidePath);
        end

        % Menu selected function: DetailsMenu
        function DetailsMenuSelected(app, ~)
            msgbox(app.DetailsInformation);
        end

        % Menu selected function: ExitMenu
        function ExitMenuSelected(app, ~)
            answer = questdlg('Are you sure you want to exit? Unsaved data will be lost.');
            if strcmp(answer,'No') || strcmp(answer, 'Cancel')
                return
            else
                close(app.UIFigure)
            end
        end

        % Menu selected function: RunMenu
        function RunMenuSelected(app, ~)
            app.RunButtonPushed();
        end

        % Value changed function: SessionSimulationsDropDown
        function SessionSimulationsDropDownValueChanged(app, ~)
            value = app.SessionSimulationsDropDown.Value;
            if strcmp(value, 'New...') == 1
                app.NumSims = app.NumSims + 1;
                app.SimNames{app.NumSims} = 'MyNewSim';
                app.SimNames{app.NumSims+1} = 'New...';
                app.Sim{app.NumSims} = Simulation([], 'MyNewSim', app.NumSims);
                app.CurrentSim = app.NumSims;
            else
                idx = find(strcmp(value, [app.SimNames]));
                app.CurrentSim = idx;
            end

            app.UpdatePlotsGUI;
            app.UpdateSimNameFields(value);
            app.UpdateImageGUI();
            app.UpdateParameterGUI();
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, ~)

            app.Sim{app.CurrentSim} = app.Sim{app.CurrentSim}.Run;
            app.TableData{app.CurrentSim, 1} = app.SimNames{app.CurrentSim};
            app.TableData{app.CurrentSim, 2} = app.Sim{app.CurrentSim}.Turnovers(1);
            app.TableData{app.CurrentSim, 3} = app.Sim{app.CurrentSim}.Turnovers(2);
            app.TableData{app.CurrentSim, 4} = app.Sim{app.CurrentSim}.Turnovers(3);

            app.ImageTableData = cell(numel({app.Sim{app.CurrentSim}.ImageFiles.name}),2);
            app.ImageTableData(:,1) = {app.Sim{app.CurrentSim}.ImageFiles.name}';
            app.ImageTableData(:,2) = {app.Sim{app.CurrentSim}.Images.Time}';

            app.UpdateDataTable;
            app.UpdateImageGUI;
            app.UpdatePlotsGUI;

        end

        % Value changed function: SessionSimulationsDropDown_4
        function SessionSimulationsDropDown_4ValueChanged(app, ~)
            value = app.SessionSimulationsDropDown_4.Value;
            if strcmp(value, 'New...') == 1
                app.SimNames{app.NumSims} = 'MyNewSim';
                app.NumSims = app.NumSims + 1;
                app.SimNames{app.NumSims} = 'New...';
                app.Sim{app.NumSims} = Simulation([], 'MyNewSim', app.NumSims);
                app.CurrentSim = app.NumSims;
            else
                idx = find(strcmp(value,app.SimNames));
                app.CurrentSim = idx;
            end

            app.UpdatePlotsGUI;
            app.UpdateSimNameFields(value);
            app.UpdateImageGUI();
            app.UpdateParameterGUI();
        end

        % Value changed function: SimulationNameEditField_4
        function SimulationNameEditField_4ValueChanged(app, ~)
            value = app.SimulationNameEditField_4.Value;
            app.UpdateSimNameFields(value)
        end

        % Value changed function: PlotDropDown
        function PlotDropDownValueChanged(app, ~)
            value = app.PlotDropDown.Value;
            app.CurrentPlotName = value;
            app.UpdatePlotsGUI;
        end

        % Value changed function: SessionSimulationsDropDown_3
        function SessionSimulationsDropDown_3ValueChanged(app, ~)
            value = app.SessionSimulationsDropDown_3.Value;
            if strcmp(value, 'New...') == 1
                app.SimNames{app.NumSims} = 'MyNewSim';
                app.NumSims = app.NumSims + 1;
                app.SimNames{app.NumSims} = 'New...';
                app.Sim{app.NumSims} = Simulation([], 'MyNewSim', app.NumSims);
                app.CurrentSim = app.NumSims;
            else
                idx = find(strcmp(value,app.SimNames));
                app.CurrentSim = idx;
            end

            app.UpdatePlotsGUI;
            app.UpdateSimNameFields(value);
            app.UpdateImageGUI();
            app.UpdateParameterGUI();
        end

        % Value changed function: SimulationNameEditField_3
        function SimulationNameEditField_3ValueChanged(app, ~)
            value = app.SimulationNameEditField_3.Value;
            app.UpdateSimNameFields(value);

        end

        % Button pushed function: PreviousPicture
        function PreviousPictureButtonPushed(app, ~)
            if isempty(app.Sim{app.CurrentSim}.Images)
                return
            elseif app.CurrentImageIndex - 1 <= 0
                return
            else
                app.CurrentImageIndex = app.CurrentImageIndex - 1;
            end
            app.UpdateImageGUI;
        end

        % Button pushed function: NextPicture
        function NextPictureButtonPushed(app, ~)
            if isempty(app.Sim{app.CurrentSim}.Images)
                return
            elseif app.CurrentImageIndex + 1 >= numel(app.Sim{app.CurrentSim}.Images)
                return
            else
                app.CurrentImageIndex = app.CurrentImageIndex + 1;
            end
            app.UpdateImageGUI;
        end

        % Cell selection callback: UITable
        function UITableCellSelection(app, event)
            indices = event.Indices;
            app.CurrentImageIndex = indices(1);
            app.UpdateImageGUI;
        end


        % Menu selected function: LoadSimulatorSessionmatMenu
        function LoadSimulatorSessionmatMenuSelected(app, ~)
            file = uigetfile('.mat', "MultiSelect",'off');
            if ~app.Sim{app.CurrentSim}.HasRun
                app.NumSims = 1;
                app.CurrentSim = app.NumSims;
                app.Sim{app.CurrentSim} = load(file, '-mat');
                app.Sim{app.CurrentSim} = app.Sim{app.CurrentSim}.this;
                app.SimNames{app.CurrentSim} = app.Sim{app.CurrentSim}.Name;
                app.SimNames{app.CurrentSim + 1} = 'New...';
            else
                app.NumSims = app.NumSims +1;
                app.CurrentSim = app.NumSims;
                app.Sim{app.CurrentSim} = load(file, '-mat', 'this');
                app.Sim{app.CurrentSim} = app.Sim{app.CurrentSim}.this;
                app.SimNames{app.CurrentSim} = app.Sim{app.CurrentSim}.Name;
                app.SimNames{app.CurrentSim + 1} = 'New...';
            end

            app.TableData{app.CurrentSim, 1} = app.SimNames{app.CurrentSim};
            app.TableData{app.CurrentSim, 2} = app.Sim{app.CurrentSim}.Turnovers(1);
            app.TableData{app.CurrentSim, 3} = app.Sim{app.CurrentSim}.Turnovers(2);
            app.TableData{app.CurrentSim, 4} = app.Sim{app.CurrentSim}.Turnovers(3);

            app.ImageTableData = cell(numel({app.Sim{app.CurrentSim}.ImageFiles.name}),2);
            app.ImageTableData(:,1) = {app.Sim{app.CurrentSim}.ImageFiles.name}';
            app.ImageTableData(:,2) = {app.Sim{app.CurrentSim}.Images.Time}';

            app.UpdateSimNameFields(app.Sim{app.CurrentSim}.Name);

            app.UpdateDataTable;
            app.UpdateImageGUI;
            app.UpdatePlotsGUI;
            app.UpdateParameterGUI;
            app.UpdateSimNameDropDowns;
        end

        % Menu selected function: LoadParametersFileMenu
        function LoadParametersFileMenuSelected(app, ~)
            [fileName, pathName] = uigetfile({'*.par', '*.txt'}, 'Select the parameter file (.par, .txt)', 'MultiSelect', 'on');
            if iscell(fileName)
                answer = questdlg('Loading multiple parameter files will overwrite existing session. Continue?', 'Warning!');
                if strcmp(answer, 'Yes') ~= 1
                    return
                end
                app.Init();
                sim = {};
                simNames = {};
                app.NumSims = numel(fileName);
                parfor i = 1:numel(fileName)
                    thoseParams = readcell(fullfile(pathName, fileName{i}), 'FileType','text', 'Delimiter', ' ');
                    sim{i} = Simulation(thoseParams(:,2), append('Sim_', num2str(i)), i);
                    simNames{i} = append('Sim_', num2str(i));

                end
                app.Sim = sim;
                app.SimNames = simNames;
                app.SimNames{app.NumSims+1} = 'New...';
            else
                app.NumSims = app.NumSims + 1;
                app.CurrentSim = app.NumSims;

                thoseParams = readcell(fullfile(pathName, fileName), 'FileType','text', 'Delimiter', ' ');
                app.Sim{app.CurrentSim} = Simulation(thoseParams(:,2), append('Sim_', num2str(app.CurrentSim)), app.CurrentSim);
                app.SimNames{app.CurrentSim} = append('Sim_', num2str(app.CurrentSim));
                app.SimNames{app.CurrentSim+1} = 'New...';
            end
            app.UpdateParameterGUI;
            app.UpdateSimNameDropDowns;

        end

        % Menu selected function: LoadDataFromFolderMenu
        function LoadDataFromFolderMenuSelected(app, ~)
            fileDir = uigetdir;

            app.Sim{app.CurrentSim}.ImageFiles = {};

            file = append(fileDir, '/Images/Image_1.jpg');
            i = 1;
            while isfile(file)

                app.Sim{app.CurrentSim}.Images(i).RGBArray = imread(file);
                app.Sim{app.CurrentSim}.Images(i).Time = int64((i-1)*2);
                app.Sim{app.CurrentSim}.ImageFiles{i} = append('Image_', num2str(i));
                i = i + 1;
                file = append(fileDir, '/Images/Image_', num2str(i), '.jpg');
            end

            app.Sim{app.CurrentSim}.Kymograph = imread(append(fileDir, '/Kymograph.jpeg'));
            app.Sim{app.CurrentSim}.InitTable  = readtable(append(fileDir, '/InitTable.csv'));
            app.Sim{app.CurrentSim}.LuminosityTable = readtable(append(fileDir, '/LuminosityTable.csv'));
            app.Sim{app.CurrentSim}.Params = readcell(append(fileDir, '/Params.csv'));
            Info = readcell(append(fileDir, '/Info.csv'), "FileType", "text", "Delimiter",',');
            app.Sim{app.CurrentSim}.Name = Info{1};
            app.Sim{app.CurrentSim}.SerialNumber = Info{2};
            app.Sim{app.CurrentSim}.Turnovers = [Info{3}, Info{4}, Info{5}];

            app.TableData{app.CurrentSim, 1} = app.SimNames{app.CurrentSim};
            app.TableData{app.CurrentSim, 2} = app.Sim{app.CurrentSim}.Turnovers(1);
            app.TableData{app.CurrentSim, 3} = app.Sim{app.CurrentSim}.Turnovers(2);
            app.TableData{app.CurrentSim, 4} = app.Sim{app.CurrentSim}.Turnovers(3);

            app.ImageTableData = cell(numel(app.Sim{app.CurrentSim}.ImageFiles),2);
            app.ImageTableData(:,1) = app.Sim{app.CurrentSim}.ImageFiles';
            app.ImageTableData(:,2) = {app.Sim{app.CurrentSim}.Images.Time}';

            app.UpdateDataTable;
            app.UpdateImageGUI;
            app.UpdateParameterGUI;
            app.UpdateSimNameDropDowns;

        end


        % Menu selected function: ExportSimulatorSessionXXXMenu
        function ExportSimulatorSessionmatMenuSelected(app, ~)
            selpath = uigetdir;
            Sims = app.Sim;
            for i = 1:app.NumSims
                Sims{i}.ExportMatFile(selpath)
            end
        end

        % Menu selected function: ExportParameterFilesMenu
        function ExportParameterFilesMenuSelected(app, ~)
            selpath = uigetdir;
            Sims = app.Sim;
            for i = 1:app.NumSims
                Sims{i}.ExportParameters(selpath);
            end
        end

        % Menu selected function: ExportRawDatatoFolderMenu
        function ExportRawDatatoFolderMenuSelected(app, ~)
            selpath = uigetdir;
            Sims = app.Sim;
            for i = 1:app.NumSims
                Sims{i}.ExportDataPackage(selpath);
            end
        end

        % Menu selected function: ExportAllPlotsMenu
        function ExportAllPlotsMenuSelected(app, ~)
            selpath = uigetdir;
            app.Sim{app.CurrentSim}.ExportAllPlots(selpath)
        end

        % Menu selected function: ExportCurrentPlotMenu
        function ExportCurrentPlotMenuSelected(app, ~)
            selpath = uigetdir;
            app.Sim{app.CurrentSim}.ExportPlot(selpath, app.CurrentPlotName);
        end


        % Drop down opening function: InitTypeDropDown
        function InitTypeDropDownOpening(app, ~)
            app.InitTypeDropDown.Value = 'RND';
            app.InitTypeDropDown.Items = {'RND', 'FIXED'};
        end

        function MaxThreadsMenuSelected(app, ~)
            answer = inputdlg(append('Threads allow for parallel computation for some operations, which increase their speed greatly. The current number of threads is: ', num2str(app.MaxThreads)), 'Set Number of Threads');
            app.MaxThreads = str2double(answer{:});
            delete(app.ParPool);
            app.ParPool = parpool(app.MaxThreads);
        end

        function CurrentMenuSelected(app, ~)
            msgbox(append('The current temporary path is: ', app.Sim{app.CurrentSim}.TempFolder));
        end

    end

    % Component initialization
    methods (Access = private)

        function app = Init(app)
            unix('rm -r ./temp')
            app.MaxThreads = 4;
            if isempty(gcp("nocreate"))
                app.ParPool = parpool(app.MaxThreads);
            end
            app.NumSims = 1;
            app.CurrentImageIndex = 1;
            app.CurrentDataTable = 1;
            app.Sim{1} = Simulation([], '', 1);
            app.CurrentSim = 1;
            app.ParameterTable.Data = app.Sim{app.CurrentSim}.Params(25:end, 2:4);
            app.ThicknessEditField.Value = app.Sim{app.CurrentSim}.Params{4, 4};
            app.WidthEditField.Value = app.Sim{app.CurrentSim}.Params{5, 4};
            app.ChangeinXEditField.Value = app.Sim{app.CurrentSim}.Params{6, 4};
            app.NumberofMTEditField.Value = app.Sim{app.CurrentSim}.Params{7, 4};
            app.NumberofXPixelsEditField.Value = app.Sim{app.CurrentSim}.Params{8, 4};
            app.NumberofYPixelsEditField.Value = app.Sim{app.CurrentSim}.Params{9, 4};
            app.MinimumTimeEditField.Value = app.Sim{app.CurrentSim}.Params{10, 4};
            app.MaximumTimeEditField.Value = app.Sim{app.CurrentSim}.Params{11, 4};
            app.FRAPTimeEditField.Value = app.Sim{app.CurrentSim}.Params{12, 4};
            app.FRAPX1EditField.Value = app.Sim{app.CurrentSim}.Params{13, 4};
            app.FRAPX2EditField.Value = app.Sim{app.CurrentSim}.Params{14, 4};
            app.FRAPY1EditField.Value = app.Sim{app.CurrentSim}.Params{15, 4};
            app.FRAPY2EditField.Value = app.Sim{app.CurrentSim}.Params{16, 4};
            app.NumberofBinsEditField.Value = app.Sim{app.CurrentSim}.Params{17, 4};
            app.OutputBMAsButton.Value = 1;
            app.PlotDropDown.Items = app.Sim{app.CurrentSim}.PlotNames;
            app.SimNames{1} = 'MyNewSim';
            app.SimNames{end+1} = 'New...';
            app.SessionSimulationsDropDown_3.Items = app.SimNames;
            app.SessionSimulationsDropDown_3.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown_2.Items = app.SimNames;
            app.SessionSimulationsDropDown_2.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown_4.Items = app.SimNames;
            app.SessionSimulationsDropDown_4.Value = app.SimNames(app.CurrentSim);
            app.SessionSimulationsDropDown.Items = app.SimNames;
            app.SessionSimulationsDropDown.Value = app.SimNames(app.CurrentSim);

            app.CurrentPlotName = app.Sim{app.CurrentSim}.PlotNames{1};

            app.InitTypeDropDown.Items = {'RND', 'FIXED'};
            app.InitTypeDropDown.Value = 'RND';

        end
        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 740 576];
            app.UIFigure.Name = 'Phragmoplast Simulator';

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.FileMenu);
            app.HelpMenu.MenuSelectedFcn = createCallbackFcn(app, @HelpMenuSelected, true);
            app.HelpMenu.Text = 'Help';

            % Create DetailsMenu
            app.DetailsMenu = uimenu(app.FileMenu);
            app.DetailsMenu.MenuSelectedFcn = createCallbackFcn(app, @DetailsMenuSelected, true);
            app.DetailsMenu.Text = 'Details';

            % Create ExitMenu
            app.ExitMenu = uimenu(app.FileMenu);
            app.ExitMenu.MenuSelectedFcn = createCallbackFcn(app, @ExitMenuSelected, true);
            app.ExitMenu.Text = 'Exit';

            % Create LoadMenu
            app.LoadMenu = uimenu(app.UIFigure);
            %             app.LoadMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadMenuSelected, true);
            app.LoadMenu.Text = 'Load';

            % Create LoadSimulatorSessionmatMenu
            app.LoadSimulatorSessionmatMenu = uimenu(app.LoadMenu);
            app.LoadSimulatorSessionmatMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadSimulatorSessionmatMenuSelected, true);
            app.LoadSimulatorSessionmatMenu.Text = 'Load Simulator Session (.mat)';

            % Create LoadParametersFileMenu
            app.LoadParametersFileMenu = uimenu(app.LoadMenu);
            app.LoadParametersFileMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadParametersFileMenuSelected, true);
            app.LoadParametersFileMenu.Text = 'Load Parameter(s) File';

            % Create LoadDataFromFolderMenu
            app.LoadDataFromFolderMenu = uimenu(app.LoadMenu);
            app.LoadDataFromFolderMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadDataFromFolderMenuSelected, true);
            app.LoadDataFromFolderMenu.Text = 'Load Data From Folder';


            % Create ExportMenu
            app.ExportMenu = uimenu(app.UIFigure);
            app.ExportMenu.Text = 'Export';



            % Create ExportParameterFilesMenu
            app.ExportParameterFilesMenu = uimenu(app.ExportMenu);
            app.ExportParameterFilesMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportParameterFilesMenuSelected, true);
            app.ExportParameterFilesMenu.Text = 'Export Parameter File(s)';

            % Create ExportCurrentPlotMenu
            app.ExportCurrentPlotMenu = uimenu(app.ExportMenu);
            app.ExportCurrentPlotMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportCurrentPlotMenuSelected, true);
            app.ExportCurrentPlotMenu.Text = 'Export Current Plot';

            % Create ExportAllPlotsMenu
            app.ExportAllPlotsMenu = uimenu(app.ExportMenu);
            app.ExportAllPlotsMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportAllPlotsMenuSelected, true);
            app.ExportAllPlotsMenu.Text = 'Export All Plots';

            % Create ExportRawDatatoFolderMenu
            app.ExportRawDatatoFolderMenu = uimenu(app.ExportMenu);
            app.ExportRawDatatoFolderMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportRawDatatoFolderMenuSelected, true);
            app.ExportRawDatatoFolderMenu.Text = 'Export Raw Data to Folder';


            app.ExportSimulatorSessionmatMenu = uimenu(app.ExportMenu);
            app.ExportSimulatorSessionmatMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportSimulatorSessionmatMenuSelected, true);
            app.ExportSimulatorSessionmatMenu.Text = 'Export Session File (.mat)';


            % Create SimulationMenu
            app.SimulationMenu = uimenu(app.UIFigure);
            app.SimulationMenu.Text = 'Simulation';

            % Create RunMenu
            app.RunMenu = uimenu(app.SimulationMenu);
            app.RunMenu.MenuSelectedFcn = createCallbackFcn(app, @RunMenuSelected, true);
            app.RunMenu.Text = 'Run';

            % Create RunSettingsMenu
            app.RunSettingsMenu = uimenu(app.SimulationMenu);
            app.RunSettingsMenu.Text = 'Run Settings';

            % Create MaxThreadsMenu
            app.MaxThreadsMenu = uimenu(app.RunSettingsMenu);
            app.MaxThreadsMenu.MenuSelectedFcn = createCallbackFcn(app, @MaxThreadsMenuSelected, true);
            app.MaxThreadsMenu.Text = 'Max Threads';

            % Create TemporaryFileLocationMenu
            app.TemporaryFileLocationMenu = uimenu(app.SimulationMenu);
            app.TemporaryFileLocationMenu.Text = 'Temporary File Location';

            % Create CurrentMenu
            app.CurrentMenu = uimenu(app.TemporaryFileLocationMenu);
            app.CurrentMenu.MenuSelectedFcn = createCallbackFcn(app, @CurrentMenuSelected, true);
            app.CurrentMenu.Text = 'Current';

            % Create ChangeRandomSeedMenu
            app.ChangeRandomSeedMenu = uimenu(app.SimulationMenu);
            app.ChangeRandomSeedMenu.MenuSelectedFcn = createCallbackFcn(app, @ChangeRandomSeedMenuSelected, true);
            app.ChangeRandomSeedMenu.Text = 'Change Random Seed';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 1 740 576];

            % Create ParametersTab
            app.ParametersTab = uitab(app.TabGroup);
            app.ParametersTab.Title = 'Parameters';

            % Create NoMatureModeButton
            app.NoMatureModeButton = uibutton(app.ParametersTab, 'state');
            app.NoMatureModeButton.ValueChangedFcn = createCallbackFcn(app, @NoMatureModeButtonValueChanged, true);
            app.NoMatureModeButton.Text = 'No Mature Mode';
            app.NoMatureModeButton.Position = [544 437 196 22];

            % Create OutputJPEGsButton
            app.OutputJPEGsButton = uibutton(app.ParametersTab, 'state');
            app.OutputJPEGsButton.ValueChangedFcn = createCallbackFcn(app, @OutputJPEGsButtonValueChanged, true);
            app.OutputJPEGsButton.Text = 'Output JPEGs';
            app.OutputJPEGsButton.Position = [544 79 89 22];

            % Create OutputBMAsButton
            app.OutputBMAsButton = uibutton(app.ParametersTab, 'state');
            app.OutputBMAsButton.ValueChangedFcn = createCallbackFcn(app, @OutputBMAsButtonValueChanged, true);
            app.OutputBMAsButton.Text = 'Output BMAs';
            app.OutputBMAsButton.Position = [649 79 91 22];



            % Create ParameterTable
            app.ParameterTable = uitable(app.ParametersTab);
            app.ParameterTable.ColumnName = {'Parameter'; 'Information'; 'Value'};
            app.ParameterTable.RowName = {};
            app.ParameterTable.ColumnSortable = [true false true];
            app.ParameterTable.ColumnEditable = [false false true];
            app.ParameterTable.CellEditCallback = createCallbackFcn(app, @ParameterTableCellEdit, true);
            app.ParameterTable.Position = [19 15 302 526];

            % Create RunButton
            app.RunButton = uibutton(app.ParametersTab, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Position = [656 481 83 22];
            app.RunButton.Text = 'Run';

            % Create NoTimeonImagesButton
            app.NoTimeonImagesButton = uibutton(app.ParametersTab, 'state');
            app.NoTimeonImagesButton.ValueChangedFcn = createCallbackFcn(app, @NoTimeonImagesButtonValueChanged, true);
            app.NoTimeonImagesButton.Text = 'No Time on Images';
            app.NoTimeonImagesButton.WordWrap = 'on';
            app.NoTimeonImagesButton.Position = [544 401 196 21];

            % Create SimulationNameEditField
            app.SimulationNameEditField = uieditfield(app.ParametersTab, 'text');
            app.SimulationNameEditField.ValueChangedFcn = createCallbackFcn(app, @SimulationNameEditFieldValueChanged, true);
            app.SimulationNameEditField.Placeholder = 'MyNewSim';
            app.SimulationNameEditField.Position = [443 480 200 22];

            % Create SimulationNameEditFieldLabel
            app.SimulationNameEditFieldLabel = uilabel(app.ParametersTab);
            app.SimulationNameEditFieldLabel.Position = [331 480 97 22];
            app.SimulationNameEditFieldLabel.Text = 'Simulation Name';

            % Create DistributionTypeDropDown
            app.DistributionTypeDropDown = uidropdown(app.ParametersTab);
            app.DistributionTypeDropDown.Items = {'Cauchy', 'Normal', 'Weibull', 'Exponential'};
            app.DistributionTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @DistributionTypeDropDownValueChanged, true);
            app.DistributionTypeDropDown.Position = [428 438 99 22];
            app.DistributionTypeDropDown.Value = 'Cauchy';

            % Create DistributionTypeDropDownLabel
            app.DistributionTypeDropDownLabel = uilabel(app.ParametersTab);
            app.DistributionTypeDropDownLabel.Position = [331 438 95 22];
            app.DistributionTypeDropDownLabel.Text = 'Distribution Type';

            % Create ThicknessEditField
            app.ThicknessEditField = uieditfield(app.ParametersTab, 'numeric');
            app.ThicknessEditField.ValueChangedFcn = createCallbackFcn(app, @ThicknessEditFieldValueChanged, true);
            app.ThicknessEditField.Position = [407 361 120 22];

            % Create ThicknessEditFieldLabel
            app.ThicknessEditFieldLabel = uilabel(app.ParametersTab);
            app.ThicknessEditFieldLabel.Position = [332 361 60 22];
            app.ThicknessEditFieldLabel.Text = 'Thickness';

            % Create WidthEditField
            app.WidthEditField = uieditfield(app.ParametersTab, 'numeric');
            app.WidthEditField.ValueChangedFcn = createCallbackFcn(app, @WidthEditFieldValueChanged, true);
            app.WidthEditField.Position = [597 361 145 22];

            % Create WidthEditFieldLabel
            app.WidthEditFieldLabel = uilabel(app.ParametersTab);
            app.WidthEditFieldLabel.HorizontalAlignment = 'right';
            app.WidthEditFieldLabel.Position = [545 361 37 22];
            app.WidthEditFieldLabel.Text = 'Width';

            % Create ChangeinXEditField
            app.ChangeinXEditField = uieditfield(app.ParametersTab, 'numeric');
            app.ChangeinXEditField.ValueChangedFcn = createCallbackFcn(app, @ChangeinXEditFieldValueChanged, true);
            app.ChangeinXEditField.Position = [418 285 110 22];

            % Create ChangeinXEditFieldLabel
            app.ChangeinXEditFieldLabel = uilabel(app.ParametersTab);
            app.ChangeinXEditFieldLabel.Position = [333 285 71 22];
            app.ChangeinXEditFieldLabel.Text = 'Change in X';

            % Create NumberofMTEditField
            app.NumberofMTEditField = uieditfield(app.ParametersTab, 'numeric');
            app.NumberofMTEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofMTEditFieldValueChanged, true);
            app.NumberofMTEditField.Position = [643 285 99 22];

            % Create NumberofMTEditFieldLabel
            app.NumberofMTEditFieldLabel = uilabel(app.ParametersTab);
            app.NumberofMTEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofMTEditFieldLabel.Position = [545 285 83 22];
            app.NumberofMTEditFieldLabel.Text = 'Number of MT';

            % Create NumberofXPixelsEditField
            app.NumberofXPixelsEditField = uieditfield(app.ParametersTab, 'numeric');
            app.NumberofXPixelsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofXPixelsEditFieldValueChanged, true);
            app.NumberofXPixelsEditField.Position = [455 324 72 22];

            % Create NumberofXPixelsEditFieldLabel
            app.NumberofXPixelsEditFieldLabel = uilabel(app.ParametersTab);
            app.NumberofXPixelsEditFieldLabel.Position = [332 324 108 22];
            app.NumberofXPixelsEditFieldLabel.Text = 'Number of X Pixels';

            % Create NumberofYPixelsEditField
            app.NumberofYPixelsEditField = uieditfield(app.ParametersTab, 'numeric');
            app.NumberofYPixelsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofYPixelsEditFieldValueChanged, true);
            app.NumberofYPixelsEditField.Position = [669 324 73 22];

            % Create NumberofYPixelsEditFieldLabel
            app.NumberofYPixelsEditFieldLabel = uilabel(app.ParametersTab);
            app.NumberofYPixelsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofYPixelsEditFieldLabel.Position = [545 324 109 22];
            app.NumberofYPixelsEditFieldLabel.Text = 'Number of Y Pixels';

            % Create MinimumTimeEditField
            app.MinimumTimeEditField = uieditfield(app.ParametersTab, 'numeric');
            app.MinimumTimeEditField.ValueChangedFcn = createCallbackFcn(app, @MinimumTimeEditFieldValueChanged, true);
            app.MinimumTimeEditField.Position = [432 243 96 22];

            % Create MinimumTimeEditFieldLabel
            app.MinimumTimeEditFieldLabel = uilabel(app.ParametersTab);
            app.MinimumTimeEditFieldLabel.Position = [332 243 85 22];
            app.MinimumTimeEditFieldLabel.Text = 'Minimum Time';

            % Create MaximumTimeEditField
            app.MaximumTimeEditField = uieditfield(app.ParametersTab, 'numeric');
            app.MaximumTimeEditField.ValueChangedFcn = createCallbackFcn(app, @MaximumTimeEditFieldValueChanged, true);
            app.MaximumTimeEditField.Position = [648 243 92 22];

            % Create MaximumTimeEditFieldLabel
            app.MaximumTimeEditFieldLabel = uilabel(app.ParametersTab);
            app.MaximumTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumTimeEditFieldLabel.Position = [545 243 88 22];
            app.MaximumTimeEditFieldLabel.Text = 'Maximum Time';

            % Create FRAPTimeEditField
            app.FRAPTimeEditField = uieditfield(app.ParametersTab, 'numeric');
            app.FRAPTimeEditField.ValueChangedFcn = createCallbackFcn(app, @FRAPTimeEditFieldValueChanged, true);
            app.FRAPTimeEditField.Position = [414 202 113 22];

            % Create FRAPTimeEditFieldLabel
            app.FRAPTimeEditFieldLabel = uilabel(app.ParametersTab);
            app.FRAPTimeEditFieldLabel.Position = [333 202 66 22];
            app.FRAPTimeEditFieldLabel.Text = 'FRAP Time';

            % Create FRAPX1EditField
            app.FRAPX1EditField = uieditfield(app.ParametersTab, 'numeric');
            app.FRAPX1EditField.ValueChangedFcn = createCallbackFcn(app, @FRAPX1EditFieldValueChanged, true);
            app.FRAPX1EditField.Position = [405 159 123 22];

            % Create FRAPX1EditFieldLabel
            app.FRAPX1EditFieldLabel = uilabel(app.ParametersTab);
            app.FRAPX1EditFieldLabel.Position = [333 159 57 22];
            app.FRAPX1EditFieldLabel.Text = 'FRAP X1 ';

            % Create FRAPY1EditField
            app.FRAPY1EditField = uieditfield(app.ParametersTab, 'numeric');
            app.FRAPY1EditField.ValueChangedFcn = createCallbackFcn(app, @FRAPY1EditFieldValueChanged, true);
            app.FRAPY1EditField.Position = [614 159 126 22];

            % Create FRAPY1EditFieldLabel
            app.FRAPY1EditFieldLabel = uilabel(app.ParametersTab);
            app.FRAPY1EditFieldLabel.HorizontalAlignment = 'right';
            app.FRAPY1EditFieldLabel.Position = [545 159 54 22];
            app.FRAPY1EditFieldLabel.Text = 'FRAP Y1';

            % Create FRAPX2EditField
            app.FRAPX2EditField = uieditfield(app.ParametersTab, 'numeric');
            app.FRAPX2EditField.ValueChangedFcn = createCallbackFcn(app, @FRAPX2EditFieldValueChanged, true);
            app.FRAPX2EditField.Position = [402 119 126 22];

            % Create FRAPX2EditFieldLabel
            app.FRAPX2EditFieldLabel = uilabel(app.ParametersTab);
            app.FRAPX2EditFieldLabel.Position = [334 119 53 22];
            app.FRAPX2EditFieldLabel.Text = 'FRAP X2';

            % Create FRAPY2EditField
            app.FRAPY2EditField = uieditfield(app.ParametersTab, 'numeric');
            app.FRAPY2EditField.ValueChangedFcn = createCallbackFcn(app, @FRAPY2EditFieldValueChanged, true);
            app.FRAPY2EditField.Position = [614 119 125 22];

            % Create FRAPY2EditFieldLabel
            app.FRAPY2EditFieldLabel = uilabel(app.ParametersTab);
            app.FRAPY2EditFieldLabel.HorizontalAlignment = 'right';
            app.FRAPY2EditFieldLabel.Position = [545 119 54 22];
            app.FRAPY2EditFieldLabel.Text = 'FRAP Y2';

            % Create NumberofBinsEditField
            app.NumberofBinsEditField = uieditfield(app.ParametersTab, 'numeric');
            app.NumberofBinsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofBinsEditFieldValueChanged, true);
            app.NumberofBinsEditField.Position = [649 202 91 22];

            % Create NumberofBinsEditFieldLabel
            app.NumberofBinsEditFieldLabel = uilabel(app.ParametersTab);
            app.NumberofBinsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofBinsEditFieldLabel.Position = [545 202 89 22];
            app.NumberofBinsEditFieldLabel.Text = 'Number of Bins';

            % Create SessionSimulationsDropDownLabel
            app.SessionSimulationsDropDownLabel = uilabel(app.ParametersTab);
            app.SessionSimulationsDropDownLabel.HorizontalAlignment = 'right';
            app.SessionSimulationsDropDownLabel.Enable = 'on';
            app.SessionSimulationsDropDownLabel.Position = [331 515 117 22];
            app.SessionSimulationsDropDownLabel.Text = 'Session Simulations:';

            % Create SessionSimulationsDropDown
            app.SessionSimulationsDropDown = uidropdown(app.ParametersTab);
            app.SessionSimulationsDropDown.ValueChangedFcn = createCallbackFcn(app, @SessionSimulationsDropDownValueChanged, true);
            app.SessionSimulationsDropDown.Enable = 'on';
            app.SessionSimulationsDropDown.Position = [463 515 257 22];

            % Create ImageBlurRadiusEditField
            app.ImageBlurRadiusEditField = uieditfield(app.ParametersTab, 'numeric');
            app.ImageBlurRadiusEditField.ValueChangedFcn = createCallbackFcn(app, @ImageBlurRadiusEditFieldValueChanged, true);
            app.ImageBlurRadiusEditField.Position = [451 79 76 22];

            % Create ImageBlurRadiusEditFieldLabel
            app.ImageBlurRadiusEditFieldLabel = uilabel(app.ParametersTab);
            app.ImageBlurRadiusEditFieldLabel.Position = [334 79 104 22];
            app.ImageBlurRadiusEditFieldLabel.Text = 'JPEG Blur Radius';

            % Create InitTypeDropDownLabel
            app.InitTypeDropDownLabel = uilabel(app.ParametersTab);
            app.InitTypeDropDownLabel.Position = [333 400 70 22];
            app.InitTypeDropDownLabel.Text = 'Init Type';

            % Create InitTypeDropDown
            app.InitTypeDropDown = uidropdown(app.ParametersTab);
            app.InitTypeDropDown.Items = {'Random', 'Fixed'};
            app.InitTypeDropDown.DropDownOpeningFcn = createCallbackFcn(app, @InitTypeDropDownOpening, true);
            app.InitTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @InitTypeDropDownValueChanged, true);
            app.InitTypeDropDown.Position = [390 400 137 22];
            app.InitTypeDropDown.Value = 'Random';

            % Create ImagesTab
            app.ImagesTab = uitab(app.TabGroup);
            app.ImagesTab.Title = 'Images';

            % Create Image
            app.Image = uiimage(app.ImagesTab);
            app.Image.Position = [349 78 355 351];

            % Create UITable
            app.UITable = uitable(app.ImagesTab);
            app.UITable.ColumnName = {'Image Number'; 'Time'};
            app.UITable.RowName = {};
            app.UITable.ColumnSortable = true;
            app.UITable.CellSelectionCallback = createCallbackFcn(app, @UITableCellSelection, true);
            app.UITable.Position = [19 15 302 524];

            % Create SimulationNameEditField_2Label
            app.SimulationNameEditField_2Label = uilabel(app.ImagesTab);
            app.SimulationNameEditField_2Label.HorizontalAlignment = 'right';
            app.SimulationNameEditField_2Label.Position = [332 451 97 22];
            app.SimulationNameEditField_2Label.Text = 'Simulation Name';

            % Create SimulationNameEditField_2
            app.SimulationNameEditField_2 = uieditfield(app.ImagesTab, 'text');
            app.SimulationNameEditField_2.Placeholder = 'MyNewSim';
            app.SimulationNameEditField_2.Position = [444 451 277 22];

            % Create SessionSimulationsDropDown_2Label
            app.SessionSimulationsDropDown_2Label = uilabel(app.ImagesTab);
            app.SessionSimulationsDropDown_2Label.HorizontalAlignment = 'right';
            app.SessionSimulationsDropDown_2Label.Enable = 'on';
            app.SessionSimulationsDropDown_2Label.Position = [331 485 117 22];
            app.SessionSimulationsDropDown_2Label.Text = 'Session Simulations:';

            % Create SessionSimulationsDropDown_2
            app.SessionSimulationsDropDown_2 = uidropdown(app.ImagesTab);
            app.SessionSimulationsDropDown_2.Enable = 'on';
            app.SessionSimulationsDropDown_2.Position = [463 485 257 22];

            % Create NumPhotos
            app.NumPhotos = uilabel(app.ImagesTab);
            app.NumPhotos.Interpreter = 'html';
            app.NumPhotos.HorizontalAlignment = 'center';
            app.NumPhotos.Position = [463 36 127 22];
            app.NumPhotos.Text = '<html><body><p><b>1 of 300</b></p></body></html>';

            % Create PreviousPicture
            app.PreviousPicture = uibutton(app.ImagesTab, 'push');
            app.PreviousPicture.ButtonPushedFcn = createCallbackFcn(app, @PreviousPictureButtonPushed, true);
            app.PreviousPicture.Position = [378 36 78 22];
            app.PreviousPicture.Text = '<<';

            % Create NextPicture
            app.NextPicture = uibutton(app.ImagesTab, 'push');
            app.NextPicture.ButtonPushedFcn = createCallbackFcn(app, @NextPictureButtonPushed, true);
            app.NextPicture.Position = [598 36 78 22];
            app.NextPicture.Text = '>>';

            % Create PlotsTab
            app.PlotsTab = uitab(app.TabGroup);
            app.PlotsTab.Title = 'Plots';

            % Create UIAxes
            app.UIAxes = uiaxes(app.PlotsTab);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [28 15 691 414];

            % Create PlotDropDownLabel
            app.PlotDropDownLabel = uilabel(app.PlotsTab);
            app.PlotDropDownLabel.HorizontalAlignment = 'right';
            app.PlotDropDownLabel.Position = [28 451 30 22];
            app.PlotDropDownLabel.Text = 'Plot:';

            % Create PlotDropDown
            app.PlotDropDown = uidropdown(app.PlotsTab);
            app.PlotDropDown.Items = {'Turnovers', 'Kymograph', 'X', 'Y', 'Z'};
            app.PlotDropDown.ValueChangedFcn = createCallbackFcn(app, @PlotDropDownValueChanged, true);
            app.PlotDropDown.Position = [73 451 100 22];
            app.PlotDropDown.Value = 'Turnovers';

            % Create SimulationNameEditField_3Label
            app.SimulationNameEditField_3Label = uilabel(app.PlotsTab);
            app.SimulationNameEditField_3Label.HorizontalAlignment = 'right';
            app.SimulationNameEditField_3Label.Position = [332 451 97 22];
            app.SimulationNameEditField_3Label.Text = 'Simulation Name';

            % Create SimulationNameEditField_3
            app.SimulationNameEditField_3 = uieditfield(app.PlotsTab, 'text');
            app.SimulationNameEditField_3.ValueChangedFcn = createCallbackFcn(app, @SimulationNameEditField_3ValueChanged, true);
            app.SimulationNameEditField_3.Placeholder = 'MyNewSim';
            app.SimulationNameEditField_3.Position = [444 451 277 22];

            % Create SessionSimulationsDropDown_3Label
            app.SessionSimulationsDropDown_3Label = uilabel(app.PlotsTab);
            app.SessionSimulationsDropDown_3Label.HorizontalAlignment = 'right';
            app.SessionSimulationsDropDown_3Label.Enable = 'on';
            app.SessionSimulationsDropDown_3Label.Position = [331 485 117 22];
            app.SessionSimulationsDropDown_3Label.Text = 'Session Simulations:';

            % Create SessionSimulationsDropDown_3
            app.SessionSimulationsDropDown_3 = uidropdown(app.PlotsTab);
            app.SessionSimulationsDropDown_3.ValueChangedFcn = createCallbackFcn(app, @SessionSimulationsDropDown_3ValueChanged, true);
            app.SessionSimulationsDropDown_3.Enable = 'on';
            app.SessionSimulationsDropDown_3.Position = [463 485 257 22];

            % Create DataTab
            app.DataTab = uitab(app.TabGroup);
            app.DataTab.Title = 'Data';

            % Create PlotDropDownLabel
            app.DataDropDownLabel = uilabel(app.DataTab);
            app.DataDropDownLabel.HorizontalAlignment = 'right';
            app.DataDropDownLabel.Position = [28 451 30 22];
            app.DataDropDownLabel.Text = 'Plot:';

            % Create PlotDropDown
            app.DataDropDown = uidropdown(app.DataTab);
            app.DataDropDown.Items = {'Turnovers', 'Luminosity', 'Init'};
            app.DataDropDown.ValueChangedFcn = createCallbackFcn(app, @DataDropDownValueChanged, true);
            app.DataDropDown.Position = [73 451 100 22];
            app.DataDropDown.Value = 'Turnovers';

            % Create UITable2
            app.UITable2 = uitable(app.DataTab);
            app.UITable2.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.UITable2.RowName = {};
            app.UITable2.ColumnSortable = true;
            app.UITable2.Position = [19 15 699 414];

            % Create SimulationNameEditField_4Label
            app.SimulationNameEditField_4Label = uilabel(app.DataTab);
            app.SimulationNameEditField_4Label.HorizontalAlignment = 'right';
            app.SimulationNameEditField_4Label.Position = [332 445 97 22];
            app.SimulationNameEditField_4Label.Text = 'Simulation Name';

            % Create SimulationNameEditField_4
            app.SimulationNameEditField_4 = uieditfield(app.DataTab, 'text');
            app.SimulationNameEditField_4.ValueChangedFcn = createCallbackFcn(app, @SimulationNameEditField_4ValueChanged, true);
            app.SimulationNameEditField_4.Placeholder = 'MyNewSim';
            app.SimulationNameEditField_4.Position = [444 445 277 22];

            % Create SessionSimulationsDropDown_4Label
            app.SessionSimulationsDropDown_4Label = uilabel(app.DataTab);
            app.SessionSimulationsDropDown_4Label.HorizontalAlignment = 'right';
            app.SessionSimulationsDropDown_4Label.Enable = 'on';
            app.SessionSimulationsDropDown_4Label.Position = [331 485 117 22];
            app.SessionSimulationsDropDown_4Label.Text = 'Session Simulations:';

            % Create SessionSimulationsDropDown_4
            app.SessionSimulationsDropDown_4 = uidropdown(app.DataTab);
            app.SessionSimulationsDropDown_4.ValueChangedFcn = createCallbackFcn(app, @SessionSimulationsDropDown_4ValueChanged, true);
            app.SessionSimulationsDropDown_4.Enable = 'on';
            app.SessionSimulationsDropDown_4.Position = [463 485 257 22];


        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SimulatorApp_Export

            % Create UIFigure and components
            createComponents(app)
            app = Init(app);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            % Register the app with App Designer
            registerApp(app, app.UIFigure)

        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end