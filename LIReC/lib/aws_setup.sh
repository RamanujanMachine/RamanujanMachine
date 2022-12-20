sudo yum -y update
sudo yum -y groupinstall "Development Tools"
sudo yum -y install openssl-devel bzip2-devel libffi-devel postgresql-devel wget
wget https://www.python.org/ftp/python/3.8.10/Python-3.8.10.tgz
tar xvf Python-3.8.10.tgz
cd Python-3.8*/
./configure --enable-optimizations
sudo make altinstall
python3.8 -V
cd ..
#sudo rm -r -f Python-* # up to you if you want cleanup or not
#then need to run pip3.8 install -e . and then you're good to go!