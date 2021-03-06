let csvToJson = require('convert-csv-to-json');
 

var express = require('express');
var path = require('path');
var bodyParser = require('body-parser');
var request = require('request');
const app = express();

const api = require('./server/routes/api');
//const fetchExperiment = require('./server/routes/fetchExperiment');

app.use(bodyParser.json());
app.use(bodyParser.urlencoded({ extended: false }));

app.use(express.static(path.join(__dirname, 'dist')));
app.use('/api',api);

app.get('/', function (req, res) {
res.sendFile(path.join(__dirname,'dist/index.html'))
});


app.listen(80, function () {
console.log('Example listening on port 3000!');
});

