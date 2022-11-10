sudo yum -y update
sudo yum -y groupinstall "Development Tools"
sudo yum -y install openssl-devel bzip2-devel libffi-devel
sudo yum -y install wget
wget https://www.python.org/ftp/python/3.10.8/Python-3.10.8.tgz
tar xvf Python-3.10.8.tgz
cd Python-3.10*/
./configure --enable-optimizations
sudo make altinstall
python3.10 -V
cd ..
sudo rm -r -f *