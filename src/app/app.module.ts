import { BrowserModule } from '@angular/platform-browser';
import { NgModule } from '@angular/core';
import {HttpClientModule} from '@angular/common/http';
import { FormsModule } from '@angular/forms';
import {ServerService} from './server.service';
import {UserInputService} from './userinput.service';
import {GlobalData} from './globaldata.service';
import {Routes, RouterModule} from "@angular/router";

import { Ng4LoadingSpinnerModule } from 'ng4-loading-spinner';

import { AppComponent } from './app.component';
import { UserInputComponent } from './user-input/user-input.component';
import { ExperimentResultsComponent} from './experiment-results/experiment-results.component';
import { AboutComponent } from './about/about.component';
import { LoadingModule } from 'ngx-loading';


const appRoutes: Routes = [
  {path: '',component: UserInputComponent},
  {path: 'experiment', component: ExperimentResultsComponent},
  {path: 'about', component: AboutComponent}
];

@NgModule({
  declarations: [
    AppComponent,
    UserInputComponent,
    ExperimentResultsComponent,
    AboutComponent
  ],
  imports: [
    Ng4LoadingSpinnerModule.forRoot() ,
    LoadingModule,
    BrowserModule,
    HttpClientModule,
    FormsModule,
    RouterModule.forRoot(appRoutes)
  ],
  providers: [ServerService,UserInputService,GlobalData],
  bootstrap: [AppComponent]
})
export class AppModule { }
