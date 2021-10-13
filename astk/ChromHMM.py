# import feature_len as fl
# from . import utils  as ul
                                                                                                                                                                                                    

def install(path):
    import os
    import zipfile
    from urllib import request

    jar = path / "ChromHMM/ChromHMM.jar" 
    if jar.exists():
        print("ChromHMM has installed!")
        exit()

    ChromHMM_url = "http://compbio.mit.edu/ChromHMM/ChromHMM.zip"
    software = path / "ChromHMM.zip"
    with request.urlopen(ChromHMM_url) as f:
        data = f.read()
        with open(software, "wb") as h:
            h.write(data)
    with zipfile.ZipFile(software, 'r') as zip_ref:
        zip_ref.extractall(path)
    os.remove(software)
    
