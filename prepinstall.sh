tar cJf myarchive.tar.xz ./GUAP
cat install.sh > GUAP_installer.sh
cat myarchive.tar.xz >> GUAP_installer.sh
chmod 777 GUAP_installer.sh
rm myarchive.tar.xz