const express = require('express');
const router = express.Router();
const csvToJson = require('csvjson');
const fs = require('fs');
const path = require('path');
var spawn = require("child_process").spawn;
const shell = require('shelljs');




router.post('/createExperiment',  (req,res) => {   
    console.log('server printing pid' + req.body.pid);
    var pid = req.body.pid;
    var expID = "";
    var pythonProcess = spawn('python',["pidDirectory.py", pid]);
    pythonProcess.stdout.on('data', function (data){
        // do something with experiment ID (data)
         expID += data;
	//shell.exec('../../bashrun/runrmutantandvasp.sh ' + pid + ' ' + expID,{shell: '/bin/bash'});
	 console.log('in server');
         var expId = {expId : expID };
	 var data;
	 var pockets;
	 var proteinInfo 
	 try {
	 	data = fs.readFileSync(path.join('/home/completePocketsNRSet',pid + 'A' + '.csv'), {encoding : 'utf8'});
		var options = {
                delimiter : " ",
                wrap : true,
                headers :"pocket area residues"};
		pockets = csvToJson.toObject(data,options);
		console.log(pockets);
		proteinInfo = {expId, pockets};
	} catch (err) {
		proteinInfo = expId;
		console.log('surf report doesnt exist');
	}
	
         res.send(JSON.stringify(proteinInfo));
    });
});



router.post('/retrieveExperiment', (req,res) => {
 	console.log(req.body);
	var info = {expId : 'experimentid recieved by server'};
	res.send(JSON.stringify(info));
 });	






module.exports = router;
