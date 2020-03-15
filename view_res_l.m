%THIS CODE DISPLAYS THE GRID PARTITIONS USING THE OUTPUT DATA
figure;


m_fmt = 0; % machine format of data - 0 -> little endian, 1 -> big endian


if(m_fmt == 0)
       
    % load the result file
    file = fopen('res');
    final_temp = fread(file, 'double');
    fclose(file);
    
    % load the parts file
    file = fopen('g_parts_l');
    g_parts = fread(file, 'uint32');
    fclose(file);
    
    % load the global vertices id order for the output file
    file = fopen('g_global_vert_ids_l');
    g_global_vert_ids = fread(file, 'uint32');
    fclose(file);
    
    % load the nodal coordinates for plotting the data
    file = fopen('g_nodal_data_l');
    g_nodal_data = fread(file, 'double');
    fclose(file);

elseif (m_fmt == 1)   
    
    % load the result file
    file = fopen('res');
    final_temp = fread(file, 'double', 'b');
    fclose(file);
    
    % load the parts file
    file = fopen('g_parts_b');
    g_parts = fread(file, 'uint32','b');
    fclose(file);
        
    % load the global vertices id order for the output file
    file = fopen('g_global_vert_ids_b');
    g_global_vert_ids = fread(file, 'uint32', 'b');
    fclose(file);
    
    % load the nodal coordinates for plotting the data
    file = fopen('g_nodal_data_b');
    g_nodal_data = fread(file, 'double', 'b');
    fclose(file);
    
end  
     
%get the colormap
cmap = jet(1024);

N_parts = max(size(g_parts)) - 1;

n_nodes = max(size(g_global_vert_ids));

x = zeros(n_nodes, 1);
y = zeros(n_nodes, 1);
c = zeros(n_nodes, 3);
counter = 1;

min_t = min(final_temp);
max_t = max(final_temp);

for i=1:n_nodes
    
   
    % global index of the node
    idx = g_global_vert_ids(i) + 1;
    
    % sort global index wise
    x(idx) = g_nodal_data(i);
    y(idx) = g_nodal_data(i + n_nodes);
    
    % calculate the color%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = final_temp(i); 
    r = ceil(((t - min_t)/(max_t - min_t)*1023 + 1));
    c(idx,:) = cmap(r, :);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
        
 sz = 1;
 scatter(x, y, sz, c);
 
 saveas(gcf,'result.png');
