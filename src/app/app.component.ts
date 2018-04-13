import { Component } from '@angular/core';
import {Router} from '@angular/router';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import { Ng4LoadingSpinnerService } from 'ng4-loading-spinner';


@Component({
  selector: 'app-root',
  templateUrl: './app.component.html',
  styleUrls: ['./app.component.css']
})
export class AppComponent {
  //loading =false;
  constructor(private serverService: ServerService,
  private router: Router,private userInputService: UserInputService,
private spinnerService: Ng4LoadingSpinnerService ){

    //this.spinnerService.show();
    this.userInputService.pidEntered.subscribe(
       (pid : string) => {
         this.spinnerService.show();
         this.submitPidToServer(pid);
          });
   }


  title = 'Protein Cavity Explorer';
  experiments = [
    {
    pid: '32',
    name: 'david'}
  ];


  onSubmit(){
    this.serverService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});
  }

  submitPidToServer(pid : string){
    this.spinnerService.show();
    console.log("submitPidToServer " + pid);
    this.serverService.createExperiment({pid: 'lake'})
    .subscribe((response) => {
      this.spinnerService.hide();
      console.log(response);
      console.log(response['pid'])},
    (error) => {console.log(error)});

this.spinnerService.hide();
    this.router.navigate(['experiment']);
  }


  submitExperimentId(experimentId : string){
    console.log(experimentId);

  }



}
