function path=set_path
    server_path='/media/kishonylab/KishonyStorage/';
    mounted_path='/home/kishonystud/kishonyserver/';
    home_path='/Users/Aharony/KishonyServer/';
    if exist(server_path)
        path=server_path; 
    elseif exist(mounted_path)
        path=mounted_path;
    elseif exist(home_path)
        path=home_path;
    end
end
