function example2_plot(S,uniform,varying,outputDirectory)
%EXAMPLE2_PLOT Plot true physical crossing times and final aggregate masses.
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
nodes=S.displayNodes; t=S.t; x=linspace(0,1,numel(nodes));
den=max([uniform.m0;uniform.mT;varying.m0;varying.mT; ...
    reshape(uniform.mnode(:,S.middle),[],1); ...
    reshape(varying.mnode(:,S.middle),[],1); ...
    reshape(S.uniform(:,S.middle),[],1);reshape(S.varying(:,S.middle),[],1)]);
blue=[0 .447 .741]; red=[.635 .078 .184];
green=[.20 .62 .20]; orange=[.88 .28 .10]; rail=[.20 .20 .20];
if strcmp(S.name,'line'), gain=.09; boundaryGain=.12; arcHeight=.075;
else, gain=.12; boundaryGain=.12; arcHeight=.10; end

fig=figure('Color','w','Units','centimeters','Position',[2 2 20 10.5], ...
    'PaperPositionMode','auto');
hold on; box off; grid on;
for k=1:numel(nodes)
    plot([x(k) x(k)],[0 1],'Color',rail,'LineWidth',1.4);
    plot(x(k),0,'o','Color',rail,'MarkerFaceColor',rail,'MarkerSize',7);
end
drawLobe(x(1),uniform.m0/den,t,boundaryGain,blue,.35,2.6,'left');
for k=2:numel(nodes)-1
    v=nodes(k);
    % NO time offsets: every curve is drawn at its true grid coordinates.
    drawLobe(x(k),uniform.mnode(:,v)/den,t,gain,green,.25,2,'right');
    plot(x(k)+gain*S.uniform(:,v)/den,t,'--','Color',green,'LineWidth',1.3);
    drawLobe(x(k),varying.mnode(:,v)/den,t,gain,orange,.25,2,'left');
    plot(x(k)-gain*S.varying(:,v)/den,t,'--','Color',orange,'LineWidth',1.3);
end
drawLobe(x(end),uniform.mT/den,t,boundaryGain,red,.35,2.6,'right');
edgeList=zeros(0,2);
for pid=1:size(S.paths,1)
    p=S.paths(pid,:); edgeList=[edgeList; p(1:end-1)' p(2:end)']; %#ok<AGROW>
end
edgeList=unique(edgeList,'rows');
for k=1:size(edgeList,1)
    left=nodes==edgeList(k,1); right=nodes==edgeList(k,2);
    drawArc(x(left),x(right),arcHeight,[.38 .38 .38],1.2);
end
for k=1:numel(nodes)
    c=[0 0 0]; if k==1, c=blue; elseif k==numel(nodes), c=red; end
    text(x(k),-.045,S.labels{k},'Interpreter','latex','Color',c,'FontSize',20, ...
        'HorizontalAlignment','center','VerticalAlignment','top');
end
h1=plot(nan,nan,'Color',green,'LineWidth',2.2);
h2=plot(nan,nan,'Color',orange,'LineWidth',2.2);
legend([h1 h2],{'Uniform capacity','Time-varying capacity'}, ...
    'Interpreter','latex','FontSize',12,'Location','northoutside', ...
    'Orientation','horizontal','Box','off');
xlim([-boundaryGain-.025 1+boundaryGain+.025]); ylim([-.13 1.03]);
set(gca,'XTick',[],'YTick',0:.2:1,'TickLabelInterpreter','latex', ...
    'FontName','Times New Roman','FontSize',14,'LineWidth',.8);
ylabel('$t$','Interpreter','latex','FontName','Times New Roman','FontSize',18);
hold off; ax=gca; ax.XAxis.Visible='off'; axtoolbar(ax,{});
if strcmp(S.name,'line'), stem='line_graph'; else, stem='simple_graph'; end
exportgraphics(fig,fullfile(outputDirectory,[stem '.png']),'Resolution',600, ...
    'BackgroundColor','white');
exportgraphics(fig,fullfile(outputDirectory,[stem '.eps']), ...
    'ContentType','vector','BackgroundColor','white');

capFig=figure('Color','w','Units','centimeters','Position',[2 2 14 8], ...
    'PaperPositionMode','auto');
hold on; profiles=zeros(S.nt,0); names={};
for k=1:numel(S.middle)
    v=S.middle(k); label=S.labels{nodes==v};
    same=find(max(abs(profiles-S.varying(:,v)),[],1)<1e-12,1);
    if isempty(same)
        profiles(:,end+1)=S.varying(:,v); names{end+1}=label; %#ok<AGROW>
    else
        names{same}=[names{same},', ',label]; %#ok<AGROW> At most six labels.
    end
end
colours=lines(size(profiles,2));
for k=1:size(profiles,2)
    plot(t,profiles(:,k),'Color',colours(k,:),'LineWidth',1.8);
end
plot(t,S.capacity*ones(S.nt,1),'k--','LineWidth',1.8); names{end+1}='Uniform';
legend(names,'Interpreter','latex','Location','eastoutside','Box','off');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('Cell capacity','Interpreter','latex','FontSize',16);
set(gca,'FontName','Times New Roman','FontSize',13,'TickLabelInterpreter','latex');
grid on; box on; hold off; axtoolbar(gca,{});
exportgraphics(capFig,fullfile(outputDirectory,[S.name '_capacity_profiles.png']), ...
    'Resolution',300,'BackgroundColor','white');
exportgraphics(capFig,fullfile(outputDirectory,[S.name '_capacity_profiles.eps']), ...
    'ContentType','vector','BackgroundColor','white');
end

function drawLobe(x0,w,t,gain,col,alph,lw,side)
if strcmp(side,'left'), x=x0-gain*w(:); else, x=x0+gain*w(:); end
% Blend with white explicitly: EPS does not support alpha transparency.
fill([x0*ones(numel(t),1);flipud(x)],[t;flipud(t)], ...
    alph*col+(1-alph)*[1 1 1],'EdgeColor','none');
plot(x,t,'Color',col,'LineWidth',lw);
end

function drawArc(xL,xR,h,col,lw)
s=linspace(0,1,80); xC=.5*(xL+xR);
plot((1-s).^2*xL+2*(1-s).*s*xC+s.^2*xR, ...
    -2*h*(1-s).*s,'Color',col,'LineWidth',lw);
end
