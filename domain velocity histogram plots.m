DomVel=[];
N_TIME_STEPS = 1428;

for i=1:numel(Allma)
    currentma=Allma{i};
    currentv=currentma.getVelocities;
    TV = vertcat( currentv{:} ); 
    V = TV(:, 2:3);
    DomVel{i}=V;
end

TQsbVy = []
TQVx = []
TQVy = []

%%Fast velocity fitting%%
for i=1:numel(DomVel)
    DV=DomVel{i}; %select domain ex: DomVel{5}=velocities for domain 21
    [Nx,Xx] = hist(DV(:,1),(N_TIME_STEPS/2)); 
    [Ny,Xy] = hist(DV(:,2),(N_TIME_STEPS/2));
    SPACE_UNITS = 'um';
    TIME_UNITS = 's';
    domain = ['Domain' num2str(i)];
    % Generate the full title
    title_strx = [domain ' 30Pa 0.7mA Vx'];
    title_stry = [domain ' 30Pa 0.7mA Vy'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maxwellian distribution
    modelFun = @(p,x) (1./(sqrt(2*pi)*p(2))).*p(1).*exp(-0.5*(((x-p(3)))./p(2)).^2); % Gaussian
    % Fit for Vx
    M_Vxstart = [6000,300,0]; % Starting values are our best initial guess for the free parametrrs p(1), p(2), etc.
    M_Vx = nlinfit(Xx,Nx,modelFun,M_Vxstart);
    M_Vx_eval = modelFun(real(M_Vx),Xx);
    % Fit for Vy
    M_Vystart = [6000,300,0]; % Starting values are our best initial guess for the free parametrrs p(1), p(2), etc.
    M_Vy = nlinfit(Xy,Ny,modelFun,M_Vystart);
    M_Vy_eval = modelFun(real(M_Vy),Xy);
    %%%%%%%%%%%%%%%%%%%%
    %%below without normalization (with or without norm does not change fitted
    %%parameters except for p(1)
    modelFun = @(p,x)  abs(p(1)).*(1./((1+((((x-p(2))).^2)./(((p(4)-1)^(-1)).*(abs(p(3)).^(2))))).^(1/(p(4)-1))));
    % Fit for Vx
    Q_Vxstart = [6000,0,300,1.2];
    Q_Vx = nlinfit(Xx,Nx,modelFun,Q_Vxstart);
    Q_Vx_eval = modelFun(real(Q_Vx),Xx);
    TQVx{i} = Q_Vx;
    % Fit for Vy
    Q_Vystart = [6000,0,300,1.8];
    Q_Vy = nlinfit(Xy,Ny,modelFun,Q_Vystart);
    Q_Vy_eval = modelFun(real(Q_Vy),Xy);
    TQVy{i} = Q_Vy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %super and super without normalization 
    modelFun = @(p,x) (abs(p(1))./((1+((((x-p(2))).^2)./(((p(4)-1)^(-1)).*(abs(p(3)).^2)))).^(1/(p(4)-1)))) + (abs(p(1))./((1+((((x-p(5))).^2)./(((p(7)-1)^(-1)).*(abs(p(6)).^2)))).^(1/(p(7)-1))));  
    % Fit for Vy
    Qsb_Vystart = [2000,0,300,1.4,0,300,1.4];
    Qsb_Vy = nlinfit(Xy,Ny,modelFun,Qsb_Vystart);
    Qsb_Vy_eval = modelFun(real(Qsb_Vy),Xy);
    TQsbVy{i} = Qsb_Vy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf;
    % Plot the histograms of x_velocities and fit different distributions
    figure
    bar(Xx,Nx)
    box off
    hold on
    plot(Xx,M_Vx_eval,'k','LineWidth',4); % Fit Gaussian distribution to the histogram
    plot(Xx,Q_Vx_eval,'m','Linewidth',2); %Fit Q Gaussian (Best fit)
    plot(0,0,'w'); %this is to create extra legend space
    title([domain ' 36Pa 0.7mA Vx']) % Change the title to reflect the parameters of the dataset used
    xlabel([ 'Velocity (' SPACE_UNITS '/' TIME_UNITS ')' ])
    ylabel('Counts')
    xlim([-2000 2000])
    set(gca,'FontSize',18); % Set axis fornsize to 18
    plot(0,0,'w'); %this is to create extra legend space
    legend({['v_{x} v_x_0=' num2str(M_Vx(3),'%.3f')], ['Mvth=' num2str(M_Vx(2))], ['q=' num2str(Q_Vx(4)) ' qvth=' num2str(abs(Q_Vx(3)), '%.3f')],['Counts ' num2str(abs(Q_Vx(1)),'%.2f')]}, 'Location', 'northeast');
    % Save the figure with the generated title as a PNG file
    % Adjust the figure size
    set(gcf, 'Position', [100, 100, 1000, 800]); % [left, bottom, width, height]
    %saveas(gcf, [title_strx '.png']);
    %saveas(gcf, [title_strx '.fig']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %y-componenet velocities
    clf;
    figure
    bar(Xy,Ny)
    box on
    hold on
    plot(Xy,M_Vy_eval,'k','LineWidth',3.5); % 
    plot(Xy,Q_Vy_eval,'m','LineWidth',2);
    plot(Xy,Qsb_Vy_eval,'g','LineWidth',2); % 
    plot(0,0,'w'); %this is to create extra legend space
    plot(0,1,'w'); 
    xlabel([ 'Velocity (' SPACE_UNITS '/' TIME_UNITS ')' ])
    ylabel('Counts')
    xlim([-2000 2000])
    set(gca,'FontSize',18); % Set axis fornsize to 18
    title([domain ' 30Pa 0.7mA Vy']) % Change the title to reflect the parameters of the dataset used
    legend({['v_{y} v_y_0=' num2str(Q_Vy(2),'%.3f')], ['M vth=' num2str(M_Vy(2),'%.3f')], ['q=' num2str(Q_Vy(4),'%.3f') ' qvth=' num2str(Q_Vy(3),'%.3f')],['qsp1=' num2str(real(Qsb_Vy(4)),'%.3f') ' qsp2=' num2str(real(Qsb_Vy(7)),'%.3f')],['qvth1=' num2str(real(Qsb_Vy(3)),'%.3f') ' qvth2=' num2str(real(Qsb_Vy(6)),'%.3f')],['qsp Counts ' num2str(real(Qsb_Vy(1)),'%.2f')]}, 'Location', 'northeast','AutoUpdate','on');
    % Save the figure with the generated title as a PNG file
    % Adjust the figure size
    set(gcf, 'Position', [100, 100, 1200, 1000]); % [left, bottom, width, height]
    %saveas(gcf, [title_stry '.png']);
    %saveas(gcf, [title_stry '.fig']);
end

% Populate matrix_TQVx3 and matrix_TQVy3 which is for thermal velocity
for i = 1:20
    row = ceil(i / 5);
    col = mod(i - 1, 5) + 1;
    matrix_TQVx3(row, col) = TQVx{i}(3);
    matrix_TQVy3(row, col) = TQVy{i}(3);
end

% Flip matrix_TQVx and matrix_TQVy horizontally
matrix_TQVx3 = flipud(matrix_TQVx3);
matrix_TQVy3 = flipud(matrix_TQVy3);

% Populate matrix_TQVx4 and matrix_TQVy4 which is for nonextensive q
for i = 1:20
    row = ceil(i / 5);
    col = mod(i - 1, 5) + 1;
    matrix_TQVx4(row, col) = TQVx{i}(4);
    matrix_TQVy4(row, col) = TQVy{i}(4);
end

% Flip matrix_TQVx4 and matrix_TQVy4 horizontally
matrix_TQVx4 = flipud(matrix_TQVx4);
matrix_TQVy4 = flipud(matrix_TQVy4);
