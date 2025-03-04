import subprocess
import requests
from os.path import basename


class UpdateChecker:
    def __init__(self,  model_full_path) -> None:
        self._tool_current_version=""
        self._model_current_version=self.get_model_num_from_filename(model_full_path)
        self._model_prefix=self.get_model_name_prefix(model_full_path)
        self._tool_latest_github=""
        self._latest_release_on_conda=False
        self._result=""

    def get_model_num_from_filename(self, filename: str) -> int:
        model_file_name=basename(filename)
        model_name=model_file_name.split(".")[0] #file name patterns is orgnanism_v#.pkl. # is the model iteration
        model_num=model_name.split("_")[-1].replace("v","")
        if model_num.isnumeric():
            return int(model_num)
        else:
            return -1
    
    def get_model_name_prefix(self, filename: str) -> str:
        model_file_name=basename(filename)
        model_name=model_file_name.split(".")[0] #file name patterns is orgnanism_v#.pkl. # is the model iteration
        model_prefix="_".join(model_name.split("_")[0:-1] )
        return model_prefix

    def get_latest_model_version(self):
        try:
            url="https://api.github.com/repos/AntonS-bio/AmpliconTyper/contents/models/"
            with requests.get(url, timeout=10) as response:
                if "status" in response.json() and response.json()["status"]=="404":
                    self._result+=f'The specified URL {url} is not available. Check with developer.'+"\n"
                    return False
                else:
                    server_model_versions=[]
                    for element in  response.json():
                        if self.get_model_name_prefix(element["name"]) == self._model_prefix:
                            server_model_versions.append(self.get_model_num_from_filename(element["name"]))
                    if max(server_model_versions)> self._model_current_version:
                        self._result+=f'Newer model version is available: {self._model_prefix}_v{max(server_model_versions)}.pkl'+"\n"
                        self._result+="You can download it here: https://github.com/AntonS-bio/AmpliconTyper/tree/main/models"+"\n"
                        return True
                    else:
                        self._result+=f'The model file {self._model_prefix}_v{self._model_current_version} is up to date'+"\n"
            return True
        except requests.exceptions.Timeout:
            self._result+="Coudn't check for update, GitHub didn't reply."+"\n"
        except requests.exceptions.ConnectionError:
            self._result+="Coudn't connect to GitHub. You may not have access to internet"+"\n"
        except requests.exceptions.InvalidURL:
            self._result+=f'The specified URL {url} is not available. Check with developer.'+"\n"
        except Exception as e:
            self._result+="Couldn't connect to GitHub."+"\n"
            self._result+=str(e)
        return False        

    def check_github_release(self) -> bool:
        try:
            url='https://api.github.com/repos/AntonS-bio/AmpliconTyper/releases/latest'
            with requests.get(url, timeout=10) as response:
                if "status" in response.json() and response.json()["status"]=="404":
                    self._result+=f'The specified URL {url} is not available. Check with developer.'+"\n"
                    return False
                else:
                    self._tool_latest_github = response.json()["tag_name"]
            return True
        except requests.exceptions.Timeout:
            self._result+="Coudn't check for update, GitHub didn't reply."+"\n"
        except requests.exceptions.ConnectionError:
            self._result+="Coudn't connect to GitHub. You may not have access to internet"+"\n"
        except requests.exceptions.InvalidURL:
            self._result+=f'The specified URL {url} is not available. Check with developer.'+"\n"
        except Exception as e:
            self._result+="Couldn't connect to GitHub. This may be due to running AmpliconTyper a lot in short time."+"\n"
        return False

    def get_installed_tool_version(self) -> bool:
        try:
            shell_stdout = subprocess.run(f'conda list | grep amplicontyper', shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
            if shell_stdout.returncode!=0:
                self._result+="Error using conda to check the currently installed version of AmpliconTyper"+"\n"
                self._result+="Run 'conda list | grep amplicontyper'. If nothing shows, you may have installed AmpliconTyper using GitHub."
                return False
            shell_stdout = shell_stdout.stdout.decode()
            lines=shell_stdout.strip().split("\n")
            if len(lines)==1:
                values=[f for f in lines[0].split(" ") if f]
                if values[0]=="amplicontyper" and len(values)>1:
                    self._tool_current_version=values[1]
                    if self._tool_current_version!=self._running_version:
                        self._result+=f'The running AmpliconTyper version is {self._running_version}'
                        self._result+=f' does not match the Conda version {self._tool_current_version}'+"\n"
                        self._result+=f'This is not a problem if you are running the tool using "python classify.py..."'+"\n"
                    return True
                else:
                    self._result+=f'Could not deteremine determine the AmpliconTyper version using conda: {lines[0]}'+"\n"
                    self._result+="Run 'conda list | grep amplicontyper'. If nothing shows, you may have installed AmpliconTyper using GitHub."
            else:
                self._result+=f'Could not deteremine determine the AmpliconTyper version using conda: {lines}'+"\n"
                self._result+="Run 'conda list | grep amplicontyper'. If nothing shows, you may have installed AmpliconTyper using GitHub."
        except Exception as e:
            self._result+="Couldn't get current version of installed tool."+"\n"
            self._result+=str(e)
        return False

    def check_bioconda_release(self) -> bool:
        try:
            result = subprocess.run(f'conda search bioconda::amplicontyper', shell=True, executable="/bin/bash", stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
            if result.returncode != 0:
                self._result+="Couldn't not check for new version of AmpliconTyper"+"\n"
            else:
                result = result.stdout.decode().strip().split("\n")
                for line in result:
                    values=[f for f in line.split(" ") if f]
                    if values[0]=="amplicontyper":
                        if len(values)>=2 and values[1]==self._tool_latest_github:
                            if self._tool_current_version==self._tool_latest_github:
                                self._result+="AmpliconTyper is up to date"+"\n"
                            else:
                                self._result+=f'Newer version ({self._tool_latest_github}) of the AmpliconTyper is available'+"\n"
                                self._result+=f'Please update using "conda install bioconda::amplicontyper={self._tool_latest_github}"'+"\n"
                            return True
                self._result+="Could not check for amplicontyper on Conda. It's probably not a problem"+"\n"
                return False
        except Exception as e:
            self._result+="Could not check for amplicontyper on Conda. It's probably not a problem"+"\n"
            self._result+=str(e)
        return False


    async def get_result(self, running_version: str) -> str:
        self._running_version=running_version
        if not self.get_installed_tool_version():
            return self._result
        if not self.check_github_release():
            return self._result
        if not self.check_bioconda_release():
            return self._result
        if not self.get_latest_model_version():
            return self._result
        return self._result
    
    @property
    def result(self) -> str:
        return self._result


