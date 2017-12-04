#!/usr/bin/env bash
set -xe

cp settings_sample.json settings.json
sed -i 's/password//' settings.json
sed -i 's/"mysqlSchema" : "swefreq"/"mysqlSchema" : "swefreq_test"/' settings.json
sed -i 's/"mysqlPort" : 3306/"mysqlPort" : 3366/' settings.json

echo "SETTINGS"
cat settings.json
echo "/SETTINGS"

mysql -u swefreq -h 127.0.0.1 -P 3366 swefreq_test < sql/swefreq.sql
mysql -u swefreq -h 127.0.0.1 -P 3366 swefreq_test < test/data/load_dummy_data.sql

cd backend

../test/01_daemon_starts.sh

python route.py --develop 1>http_log.txt 2>&1 &
BACKEND_PID=$!

function exit_handler() {
    kill -9 $BACKEND_PID

    echo "THE HTTP LOG WAS:"
    cat http_log.txt
}

trap exit_handler EXIT

sleep 2 # Lets wait a little bit so the server has started
python test.py -v
